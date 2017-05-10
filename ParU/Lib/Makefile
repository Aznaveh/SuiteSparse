#===============================================================================
# SuiteSparseQR/Lib/Makefile
#===============================================================================

LIBRARY = libspqr
VERSION = 2.0.7
SO_VERSION = 2

default: library

ccode: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

# SPQR depends on CHOLMOD, AMD, COLAMD, LAPACK, the BLAS and SuiteSparse_config
LDLIBS += -lamd -lcolamd -lcholmod -lsuitesparseconfig $(LAPACK) $(BLAS)

# compile and install in SuiteSparse/lib
library:
	$(MAKE) install INSTALL=$(SUITESPARSE)

#-------------------------------------------------------------------------------

all: library

purge: distclean

distclean: clean
	- $(RM) -r $(PURGE)

clean:
	- $(RM) -r $(CLEAN)

INC = ../Include/spqr.hpp ../Include/SuiteSparseQR_C.h \
	../Include/SuiteSparseQR_definitions.h \
	../Include/SuiteSparseQR.hpp

OBJ = \
#    spqr_rmap.o \
#    SuiteSparseQR_C.o \
#    SuiteSparseQR_expert.o \
#    spqr_parallel.o \
#    spqr_kernel.o \
#    spqr_analyze.o \
#    spqr_assemble.o \
#    spqr_cpack.o \
#    spqr_csize.o \
#    spqr_fcsize.o \
#    spqr_debug.o \
#    spqr_front.o \
#    spqr_factorize.o \
#    spqr_freenum.o \
#    spqr_freesym.o \
#    spqr_freefac.o \
#    spqr_fsize.o \
#    spqr_maxcolnorm.o \
#    spqr_rconvert.o \
#    spqr_rcount.o \
#    spqr_rhpack.o \
#    spqr_rsolve.o \
#    spqr_stranspose1.o \
#    spqr_stranspose2.o \
#    spqr_hpinv.o \
#    spqr_1fixed.o \
#    spqr_1colamd.o \
#    SuiteSparseQR.o \
#    spqr_1factor.o \
#    spqr_cumsum.o \
#    spqr_shift.o \
#    spqr_happly.o \
#    spqr_panel.o \
#    spqr_happly_work.o \
#    SuiteSparseQR_qmult.o \
#    spqr_trapezoidal.o \
#    spqr_larftb.o \
#    spqr_append.o \
#    spqr_type.o \
#    spqr_tol.o
#
$(OBJ): $(INC)

I = -I../../CHOLMOD/Include -I../../SuiteSparse_config -I../Include

C = $(CXX) $(CF) $(SPQR_CONFIG) $(CHOLMOD_CONFIG) $(CONFIG_PARTITION) $(I)
#-------------------------------------------------------------------------------

static: $(AR_TARGET)

$(AR_TARGET): $(OBJ)
	$(ARCHIVE) $@ $^
	- $(RANLIB) $@

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

# install SPQR
install: $(AR_TARGET) $(INSTALL_LIB)/$(SO_TARGET)

$(INSTALL_LIB)/$(SO_TARGET): $(OBJ)
	@mkdir -p $(INSTALL_LIB)
	@mkdir -p $(INSTALL_INCLUDE)
	@mkdir -p $(INSTALL_DOC)
	$(CXX) $(SO_OPTS) $^ -o $@ $(LDLIBS)
	( cd $(INSTALL_LIB) ; ln -sf $(SO_TARGET) $(SO_PLAIN) )
	( cd $(INSTALL_LIB) ; ln -sf $(SO_TARGET) $(SO_MAIN) )
	$(CP) ../Include/SuiteSparseQR.hpp $(INSTALL_INCLUDE)
	$(CP) ../Include/SuiteSparseQR_C.h $(INSTALL_INCLUDE)
	$(CP) ../Include/SuiteSparseQR_definitions.h $(INSTALL_INCLUDE)
	$(CP) ../Include/spqr.hpp $(INSTALL_INCLUDE)
	$(CP) ../Doc/spqr_user_guide.pdf $(INSTALL_DOC)
	$(CP) ../README.txt $(INSTALL_DOC)/SPQR_README.txt
	chmod 755 $(INSTALL_LIB)/$(SO_TARGET)
	chmod 644 $(INSTALL_INCLUDE)/SuiteSparseQR.hpp
	chmod 644 $(INSTALL_INCLUDE)/SuiteSparseQR_C.h
	chmod 644 $(INSTALL_INCLUDE)/SuiteSparseQR_definitions.h
	chmod 644 $(INSTALL_INCLUDE)/spqr.hpp
	chmod 644 $(INSTALL_DOC)/spqr_user_guide.pdf
	chmod 644 $(INSTALL_DOC)/SPQR_README.txt

# uninstall SPQR
uninstall:
	$(RM) $(INSTALL_LIB)/$(SO_TARGET)
	$(RM) $(INSTALL_LIB)/$(SO_PLAIN)
	$(RM) $(INSTALL_LIB)/$(SO_MAIN)
	$(RM) $(INSTALL_INCLUDE)/SuiteSparseQR.hpp
	$(RM) $(INSTALL_INCLUDE)/SuiteSparseQR_C.h
	$(RM) $(INSTALL_INCLUDE)/SuiteSparseQR_definitions.h
	$(RM) $(INSTALL_INCLUDE)/spqr.hpp
	$(RM) $(INSTALL_DOC)/spqr_user_guide.pdf
	$(RM) $(INSTALL_DOC)/SPQR_README.txt

