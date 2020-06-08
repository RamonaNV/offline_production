ifdef I3_PORTS
dir	=	-I$(I3_PORTS)/include -L$(I3_PORTS)/lib
else
dir	=
endif
gsl	=	$(dir) -lgsl -lgslcblas

ifdef opcl
opts	=	-I$(opcl)/include -L$(opcl)/lib/x86_64 -lOpenCL
else
opts	=	-framework OpenCL
endif

gcpp	=	$(CXX) -O2 #--fast-math
nvcc	=	nvcc -O2 --use_fast_math

all:	ocl

ocl:
	$(gcpp) llh.cxx -o llh $(gsl) ../ocl/ppc.o $(opts)

inv:
	c++ inv.cxx -o inv $(gsl)

ico:
	c++ ico.cxx -o ico
