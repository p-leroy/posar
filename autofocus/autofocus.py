#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 15:10:30 2020

@author: pleroy
"""

# <IMPORT>

import logging, os, sys

# Deal with minor differences between PySide2 and PyQt5
try:
    from PySide2 import QtCore, QtGui, QtWidgets
    Signal = QtCore.Signal
    Slot = QtCore.Slot
except ImportError:
    from PyQt5 import QtCore, QtGui, QtWidgets
    Signal = QtCore.pyqtSignal
    Slot = QtCore.pyqtSlot

from mainwindow import *
sys.path.insert(0, "/home/pleroy/DEV/processing/PoSAR-MC")
import posarmctools.epsgtools as epsg
import posarmctools.posar as posar
import posarmctools.sbg as sbg

# </IMPORT>

epsg3xxx = epsg.epsg3948
epsgStr = "epsg3948"
expectedPeriod = 1000000

class Signaller(QtCore.QObject):
    signal = Signal(str, logging.LogRecord)

class QtHandler(logging.Handler):
    def __init__(self, slotfunc, *args, **kwargs):
        super(QtHandler, self).__init__(*args, **kwargs)
        self.signaller = Signaller()
        self.signaller.signal.connect(slotfunc)

    def emit(self, record):
        s = self.format(record)
        self.signaller.signal.emit(s, record)
        QtWidgets.QApplication.processEvents()

class MyMainWindow(QtWidgets.QMainWindow):
    
    logging_colors = {
        logging.DEBUG: QtGui.QColorConstants.Black,
        logging.INFO: QtGui.QColorConstants.Blue,
        logging.WARNING: QtGui.QColorConstants.Magenta,
        logging.ERROR: QtGui.QColorConstants.Red,
        logging.CRITICAL: QtGui.QColorConstants.Magenta,
    }
    
    def __init__(self, parent=None):
        super(MyMainWindow, self).__init__(parent)
        
        self.logs = None
        
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.textEdit = self.ui.textEdit
        
        self.logger = logging.getLogger()
        
        self.handler = QtHandler(self.append_log)
        # format string
        fs = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        formatter = logging.Formatter(fs)
        self.handler.setFormatter(formatter)
        self.logger.addHandler(self.handler)
        
        self.readSettings()
        
        self.buildView()

        self.ui.pushButton_home.setText(self.sbg_dir)
        self.ui.pushButton_home.clicked.connect(self.on_click)
        self.ui.pushButton_readLogs.clicked.connect(self.read_sbg_logs)
        self.ui.pushButton_readPoSAR.clicked.connect(self.read_posar_data)
        self.ui.pushButton_clear_logs.clicked.connect(self.clear_logs)
        self.ui.viewHours.itemDoubleClicked.connect(self.read_hour)
        
        self.ui.dockWidget_2.show()
        self.ui.dockWidget_8.show()

    def closeEvent(self, event=None):
        self.writeSettings()
        super().closeEvent(event)

    def writeSettings(self):
        settings = QtCore.QSettings("IETR", "autofocus")
        settings.setValue("MyMainWindow/geometry", self.saveGeometry())
        settings.setValue("MyMainWindow/windowState", self.saveState())
        settings.setValue("sbg_dir", self.sbg_dir)

    def readSettings(self):
        settings = QtCore.QSettings("IETR", "autofocus")
        self.restoreGeometry(settings.value("MyMainWindow/geometry", self.saveGeometry()))
        self.restoreState(settings.value("MyMainWindow/windowState", self.saveState()))
        self.sbg_dir = settings.value("sbg_dir", QtCore.QDir.homePath())
            
    def buildView(self):
        self.modelDirs = QtWidgets.QFileSystemModel()
        self.modelDirs.setRootPath(self.sbg_dir)
        
        self.treeView = self.ui.treeView
        self.treeView.setModel(self.modelDirs)
        self.treeView.setRootIndex(self.modelDirs.index(self.sbg_dir))
        self.treeView.setColumnHidden(1, True)
        self.treeView.setColumnHidden(2, True)
        self.treeView.setColumnHidden(3, True)
        
        self.viewHours = self.ui.viewHours
        
        self.treeView.doubleClicked.connect(self.treeview_doubleclicked)

    def update_sbg_dir(self):
        self.modelDirs.setRootPath(self.sbg_dir)
        self.treeView.setRootIndex(self.modelDirs.index(self.sbg_dir))
        self.ui.pushButton_home.setText(self.sbg_dir)

    @Slot(bool)
    def on_click(self, checked):
        fileDialog = QtWidgets.QFileDialog()
        fileDialog.setFileMode(QtWidgets.QFileDialog.Directory)
        dir_ = fileDialog.getExistingDirectory(self, "Home directory", self.sbg_dir)
        if dir_ != "":
            self.sbg_dir = dir_
            self.update_sbg_dir()

    @Slot(QtCore.QModelIndex)
    def treeview_doubleclicked(self, idx):
        path = self.modelDirs.filePath(idx) 
        if path[-11:] == "_flight.ini":
            self.flight = posar.FlightParams(path)
            with open(path) as f:
                self.textEdit.clear()
                for line in f.readlines():
                    self.textEdit.append(line)
                    
            for k, v in self.flight.hours.items():
                item = QtWidgets.QListWidgetItem(f"{k} {v[0]} {v[1]}", self.viewHours)
                
            filename = os.path.join(self.flight.dir_posar, self.flight.refs)
            self.pts, self.ptsEpsg, self.ptsDict = epsg.getReferencePoints(filename, epsg3xxx, file="ini")
            self.ui.pushButton_readLogs.setEnabled(True)

    @Slot(QtWidgets.QListWidgetItem)
    def read_hour(self, item):
        if self.logs is not None:
            self.ui.textEdit_logging.append("")
            text = item.text()
            hour = text.split(" ")[0]
            rec_date = f"{self.flight.day}_{hour}"
            rec_dir = os.path.join(self.flight.dir_posar, self.flight.day, rec_date)
            self.record = posar.Record(rec_dir, "record", "bin", utc=self.logs.utc, period=expectedPeriod, version="X_v1")
        else:
            self.logger.error("self.logs does not exists")

    @Slot(str, logging.LogRecord)
    def append_log(self, msg, record):
        color = self.logging_colors.get(record.levelno, QtGui.QColor.black)
        self.ui.textEdit_logging.setTextColor(color)
        self.ui.textEdit_logging.append(msg)
        
    @Slot(bool)
    def read_sbg_logs(self, checked):
        self.logger.info("read_sbg_logs")
        self.logs = sbg.Logs(self.flight, epsg3xxx)
        self.ui.pushButton_readPoSAR.setEnabled(True)
    
    @Slot(bool)
    def clear_logs(self, checked):
        self.ui.textEdit_logging.clear()
    
    @Slot(bool)
    def read_posar_data(self, checked):
        self.logger.info("read_posar_data")
        self.rec_dates = [f"{self.flight.day}_" + h for h in self.flight.hours]
        self.rec_dirs = [os.path.join(self.flight.dir_posar, self.flight.day, d) for d in self.rec_dates]
        self.records = [posar.Record(rec_dir, "record", "bin", utc=self.logs.utc, period=expectedPeriod, version="X_v1") for rec_dir in self.rec_dirs]

    
def main():
        logging.getLogger().setLevel(logging.INFO)
        if not QtWidgets.QApplication.instance():
            print("NOT QApplication.instance()")
            app = QtWidgets.QApplication(sys.argv)
        else:
            print("QApplication.instance()")
            app = QtWidgets.QApplication.instance()
        mainWindow = MyMainWindow()
        mainWindow.show()
        app.exec()

if __name__ == "__main__":    
    main()
