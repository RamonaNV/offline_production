__all__=['glam', 'spglam', 'photonics', 'splinetable']

try:
	import numpy
	from icecube.load_pybindings import load_pybindings
	load_pybindings(__name__, __path__)
	del numpy
except ImportError:
	pass

try:
	from icecube.photospline import splinetable
	from icecube.photospline import splinefitstable
	from icecube.photospline import glam
	try:
		from icecube.photospline import spglam
	except ImportError:
		pass
	from icecube.photospline import numpy_extensions
	from icecube.photospline import photonics
except ImportError:
	pass

