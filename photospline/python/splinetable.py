class SplineTable(object):
	def __init__(self):
		# Context information for the fit
		self.order = 2
		self.knots = []
		self.periods = []

		# The tensor-product basis function coefficients
		self.coefficients = None

		# logarithmic bias
		self.bias = 0

		# extent of supported region
		self.extents = []

		# parity
		self.parity = 1

		# geometry type
		self.geometry = 2
		# geometry class of photonics table (0 is point souce,
		# 1 is inf. muon)
		self.geotype = 0

		# level of photonics table
		self.level = 1

		# group phase velocity if spline is fit
		# from a photonics table
		self.ngroup = -1.0
