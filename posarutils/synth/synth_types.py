class GEOMt(object):
	"""docstring fos GEOMt"""
	def __init__(self):
		super(GEOMt, self).__init__()

		self.gr_sr = 0 # 0:ground range, 1:slant_range DO NOT CHANGE !!!!

		# GROUND range and AZIMUTH focusing domain
		self.gr_pos_min = 0
		self.gr_pos_max = 0
		self.az_pos_min = 0
		self.az_pos_max = 0
		self.el_pos = 0
		# Required SLANT range and AZIMUTH resolution
		# (automatically set to max(res,res_min) with res_min the system min resolution)
		self.rg_res     = 0
		self.az_res     = 0
		# SLANT OR GROUND range and AZIMUTH pixel spacing
		self.d_gr       = 0
		self.d_az       = 0

	def print(self):
		print( 'gr_sr ' + str(self.gr_sr ) )

		print( 'gr_pos_min ' + str(self.gr_pos_min) )
		print( 'gr_pos_max ' + str(self.gr_pos_max) )
		print( 'az_pos_min ' + str(self.az_pos_min) )
		print( 'az_pos_max ' + str(self.az_pos_max) )

		print( 'el_pos ' + str(self.el_pos) )
		print( 'rg_res ' + str(self.rg_res) )
		print( 'az_res ' + str(self.az_res) )
		print( 'd_gr ' + str(self.d_gr) )
		print( 'd_az ' + str(self.d_az) )

class CFGt(object):
	"""docstring for CFGt"""
	def __init__(self):
		super(CFGt, self).__init__()

		self.Npos = 0
		self.Nf = 0
		self.Fmin = 0
		self.Fmax = 0
		self.c = 0
		self.kc = 0
		self.d_shift = 0
		self.n_up_smp = 0
		self.up_smp = 0

class PRMt(object):
	"""docstring for PRMt"""
	def __init__(self):
		super(PRMt, self).__init__()

		self.up_smp = 0
		self.check = 0
		self.xamin = 0
		self.xamax = 0
		self.force_res = 0
		self.proc_Dphi = 0
		self.min_ant_pat = 0
		self.sin_phi_max = 0

	def print(self):
		print( 'up_smp ' + str(self.up_smp ) )
		print( 'check ' + str(self.check ) )
		print( 'xamin ' + str(self.xamin ) )
		print( 'xamax ' + str(self.xamax ) )
		print( 'force_res ' + str(self.force_res ) )
		print( 'proc_Dphi ' + str(self.proc_Dphi ) )
		print( 'min_ant_pat ' + str(self.min_ant_pat ) )
		print( 'sin_phi_max ' + str(self.sin_phi_max ) )
