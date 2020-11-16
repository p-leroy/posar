import numpy as np

class AntPos(object):
	"""docstring for AntPosT"""
	def __init__(self, linspaceStart, linspaceStop, nbPoints, offset, altitude):
		super(AntPos, self).__init__()

		self.nbPoints = nbPoints
		self.altitude = altitude
		self.offset = offset

		self.Tx_x = np.linspace( linspaceStart, linspaceStop, nbPoints ) + offset / 2
		self.Tx_y = self.Tx_x * 0
		self.Tx_z = self.Tx_x * 0 + altitude
		self.Rx_x = np.linspace( linspaceStart, linspaceStop, nbPoints ) - offset / 2
		self.Rx_y = self.Rx_x * 0
		self.Rx_z = self.Rx_x * 0 + altitude
		self.mean_x = ( self.Tx_x + self.Rx_x ) / 2
		self.mean_y = ( self.Tx_y + self.Rx_y ) / 2
		self.mean_z = ( self.Tx_z + self.Rx_z ) / 2

class AntPos3D(object):
	"""docstring for AntPosT"""
	def __init__(self, firstRamp, lastRamp, xyz, offset, altitude):
		super(AntPos3D, self).__init__()

		self.altitude = altitude
		self.offset = offset

		self.Tx_x = xyz[ firstRamp : lastRamp + 1, 2] + offset / 2
		self.Tx_y = xyz[ firstRamp : lastRamp + 1, 3]
		self.Tx_z = self.Tx_x * 0 + altitude
		self.Rx_x = xyz[ firstRamp : lastRamp + 1, 2] - offset / 2
		self.Rx_y = xyz[ firstRamp : lastRamp + 1, 3]
		self.Rx_z = self.Rx_x * 0 + altitude
		self.mean_x = xyz[ firstRamp : lastRamp + 1, 2]
		self.mean_y = xyz[ firstRamp : lastRamp + 1, 3]
		self.mean_z = ( self.Tx_z + self.Rx_z ) / 2

class AntPosGPS(object):
	"""docstring for AntPosT"""
	def __init__(self, firstRamp, lastRamp, rampNumber_timeStamp_xyz, offset):
		super(AntPosGPS, self).__init__()

		self.offset = offset
		# Tx
		self.Tx_x = rampNumber_timeStamp_xyz[ firstRamp : lastRamp, 2] + offset / 2
		self.Tx_y = rampNumber_timeStamp_xyz[ firstRamp : lastRamp, 3]
		self.Tx_z = rampNumber_timeStamp_xyz[ firstRamp : lastRamp, 4]
		# Rx
		self.Rx_x = rampNumber_timeStamp_xyz[ firstRamp : lastRamp, 2] - offset / 2
		self.Rx_y = rampNumber_timeStamp_xyz[ firstRamp : lastRamp, 3]
		self.Rx_z = rampNumber_timeStamp_xyz[ firstRamp : lastRamp, 4]
		# mean
		self.mean_x = rampNumber_timeStamp_xyz[ firstRamp : lastRamp, 2]
		self.mean_y = rampNumber_timeStamp_xyz[ firstRamp : lastRamp, 3]
		self.mean_z = rampNumber_timeStamp_xyz[ firstRamp : lastRamp, 4]

class AntPosPoSARGB(object):
	"""docstring for AntPosT"""
	def __init__(self, Txx, Txy, Txz, Rxx, Rxy, Rxz):
		super(AntPos, self).__init__()

		self.Txx = Txx
		self.Txy = Txy
		self.Txz = Txz
		self.Rxx = Rxx
		self.Rxy = Rxy
		self.Rxz = Rxz
		self.mean_x = ( self.Txx + self.Rxx ) / 2
		self.mean_y = ( self.Txy + self.Rxy ) / 2
		self.mean_z = ( self.Txz + self.Rxz ) / 2
