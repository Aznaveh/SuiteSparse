#===============================================================================
# SuiteSparseQR/Demo/Makefile
#===============================================================================

default: all

ccode: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

CLIB = $(LDFLAGS) -L../../lib -lspqr -lsuitesparseconfig -lcholmod -lamd \
        -lcolamd $(LIB_WITH_PARTITION)  $(LDLIBS)

# use the BLAS and LAPACK defined by SuiteSparse_config.mk; do not use valgrind 
FLIB = $(LAPACK) $(BLAS)
V =

# To use Valgrind and the plain BLAS and plain LAPACK (non-optimized):
# FLIB = -lgfortran -llapack_plain -lblas_plain -lg2c
# V = valgrind --quiet

all: library qrsimple 
	- $(V) ./qrsimple < ../Matrix/ash219.mtx
	- $(V) ./qrsimple < ../Matrix/west0067.mtx

library: metis
	( cd ../Lib ; $(MAKE) )
	( cd ../../AMD ; $(MAKE) library )
	( cd ../../SuiteSparse_config ; $(MAKE) library )
	- ( cd ../../CHOLMOD && $(MAKE) library )
	- ( cd ../../COLAMD && $(MAKE) library )
	- ( cd ../../CCOLAMD && $(MAKE) library )
	- ( cd ../../CAMD && $(MAKE) library )
metis: ../../include/metis.h

../../include/metis.h:
	- ( cd ../.. && $(MAKE) metis )

purge: distclean

distclean: clean
	- $(RM) qrdemo qrdemo_gpu qrdemoc qrsimple qrsimplec X.mtx
	- $(RM) R.mtx C.mtx E.txt gpu_results.txt qrdemo_gpu2 qrdemo_gpu3
	- $(RM) *.dot pfile tfile
	- $(RM) -r $(PURGE)

clean:
	- $(RM) -r $(CLEAN)

INC = ../Include/spqr.hpp ../Include/SuiteSparseQR_C.h \
	../Include/SuiteSparseQR_definitions.h \
	../Include/SuiteSparseQR.hpp Makefile

I = -I../../include 

C = $(CXX) $(CF) $(SPQR_CONFIG) $(CONFIG_PARTITION) $(I) \
	$(CHOLMOD_CONFIG)

LIBS = $(CLIB) $(FLIB) $(TBB)


qrsimple: qrsimple.cpp $(INC)
	$(C) qrsimple.cpp -o qrsimple $(LIBS)
