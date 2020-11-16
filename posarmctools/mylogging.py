import json
import os.path
import re
import ipykernel
import requests
import logging

from urllib.parse import urljoin
from notebook.notebookapp import list_running_servers

def get_notebook_name():
    """
    Return the full path of the jupyter notebook.
    """
    kernel_id = re.search('kernel-(.*).json',
        ipykernel.connect.get_connection_file()).group(1)
    servers = list_running_servers()
    for ss in servers:
        response = requests.get(urljoin(ss['url'], 'api/sessions'),
            params={'token': ss.get('token', '')})
        for nn in json.loads(response.text):
            if nn['kernel']['id'] == kernel_id:
                relative_path = nn['notebook']['path']
                return os.path.join(ss['notebook_dir'], relative_path)

def loggingSetFh(name, logger, level=logging.INFO):
    # create file handler which logs even debug messages
    fh = logging.FileHandler(name, mode='w')
    fh.setLevel(level)
    # set the formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

def getNotebookLogger(level=logging.INFO):
    modulename = get_notebook_name().split('/')[-1].split('.ipynb')[0]
    logger = logging.getLogger()
    logger.setLevel(level)
    logfile = modulename + ".log"
    loggingSetFh(logfile, logger, level)
    return logger

def htmlGenerator(name):
    logFile = f"{name}.log"
    htmlFile = f"{name}.html"
    
    body = "body { background-color:white; }\n"
    error = ".error { color:#9C0006;background-color:#FFC7CE; }\n"
    warning = ".warning { color:#9C6500; background-color:#FFEB9C; }\n"
    info = ".info { color:#006100; background-color:#C6EFCE; }\n"
    record = ".record { color:#00008B; background-color:#87CEFA; }\n"
    other = ".other { background-color:white; }\n"
    style = "<style>\n" + body + error + warning + record + info + other + "</style>"
    
    with open(htmlFile, "w+") as htmlOut:
        htmlOut.write(style + "\n")
        with open(logFile, "r") as file:
            for num, line in enumerate(file.readlines()):
                if " - ERROR - " in line:
                    class_ = "error"
                elif " - WARNING - " in line:
                    class_ = "warning"
                elif " - RECORD - " in line:
                    class_ = "record"
                elif " - INFO - " in line:
                    class_ = "info"
                else:
                    class_ = "other"
                start = f"<span class=\"{class_}\">"
                stop = f"</span>"
                htmlOut.write(start + line + stop + "<br>\n")
