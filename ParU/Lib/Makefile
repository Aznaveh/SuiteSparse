#===============================================================================
# ParU/Lib/Makefile
#===============================================================================

LIBRARY = libparu
VERSION = 0.0.1
SO_VERSION = 0

default: library

ccode: all

include ../../SuiteSparse_config/SuiteSparse_config.mk


# ParU depends on SPQR, CHOLMOD, AMD, COLAMD, LAPACK, the BLAS and 
# SuiteSparse_config
#LDLIBS += -lspqr -lamd -lcolamd -lcholmod -lsuitesparseconfig $(LAPACK) $(BLAS)
#
# newer version  of ParU depends on UMFPACK, CHOLMOD, AMD, COLAMD  
# the BLAS and SuiteSparse_config
LDLIBS += -lumfpack -lamd -lcolamd -lcholmod -lsuitesparseconfig \
		  $(LAPACK) $(BLAS)

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

INC = ../Include/ParU.hpp 

C11 = -std=c++11
CC11 = -std=c11

OBJ = \
    paru_mem.o \
    paru_memset.o \
    paru_memcpy.o \
    paru_heap.o \
    paru_hash.o \
    paru_pivotal.o \
    paru_analyze.o \
	paru_init_rowFronts.o\
	paru_tuples.o\
	paru_front.o\
	paru_factorize.o\
	paru_fs_factorize.o\
	paru_create_element.o\
	paru_assemble_row2U.o\
	paru_trsm.o\
	paru_dgemm.o\
	paru_print.o\
	paru_init_rel.o\
	paru_update_rowDeg.o\
	paru_update_rel_ind.o\
	paru_write.o\
	paru_cumsum.o\
	paru_intersection.o\
	paru_bin_search.o\
	paru_assemble.o\
	paru_prior_assemble.o \
	paru_perm.o \
	paru_lsolve.o \
	paru_usolve.o \
	paru_solve.o \
	paru_residual.o \
	paru_backward.o \
	paru_gaxpy.o \
	paru_norms.o \
	paru_full_summed.o \

$(OBJ): $(INC)

I = -I ../../UMFPACK/Source -I../../UMFPACK/Include -I ../../CHOLMOD/Include\
	-I ../../AMD/Include  -I ../../SuiteSparse_config -I../Include

FLAG = -Wno-write-strings   -Wall # -Werror -Wnull-dereference -Wl,--verbose\
         -Wextra  -Wpedantic -Wshadow -Wnon-virtual-dtor #development


ifeq ($(DEBUG),yes)
	#to compile with -g options use 'make DBUG=yes'
	FLAG = -Wno-write-strings -g
endif

ifeq ($(PROF),yes)
	#to compile with -g options use 'make DBUG=yes'
	FLAG = -Wno-write-strings -no-pie -g -pg 
endif



#C = $(CXX) $(CF) $(SPQR_CONFIG) $(CHOLMOD_CONFIG) $(CONFIG_PARTITION) $(I)
C = $(CXX) $(CF)  $(CHOLMOD_CONFIG) $(CONFIG_PARTITION) $(I) $(FLAG) $(C11)
C2 = $(CC) $(CF)  $(CHOLMOD_CONFIG) $(CONFIG_PARTITION) $(I) $(FLAG) $(CC11)
#-------------------------------------------------------------------------------

static: $(AR_TARGET)

$(AR_TARGET): $(OBJ)
	$(ARCHIVE) $@ $^
	- $(RANLIB) $@

#-------------------------------------------------------------------------------

paru_mem.o: ../Source/paru_mem.cpp
	$(C) -c $<

paru_memset.o: ../Source/paru_memset.cpp
	$(C) -c $<

paru_memcpy.o: ../Source/paru_memcpy.cpp
	$(C) -c $<

paru_heap.o: ../Source/paru_heap.cpp
	$(C) -c $<

paru_hash.o: ../Source/paru_hash.cpp
	$(C) -c $<
	
paru_pivotal.o: ../Source/paru_pivotal.cpp
	$(C) -c $<

paru_analyze.o: ../Source/paru_analyze.cpp
	$(C) -c $<

paru_init_rowFronts.o: ../Source/paru_init_rowFronts.cpp
	$(C) -c $<

paru_tuples.o: ../Source/paru_tuples.cpp
	$(C) -c $<

#paru_assemble.o: ../Source/paru_assemble.cpp
#	$(C) -c $<

paru_front.o: ../Source/paru_front.cpp
	$(C) -c $<

paru_factorize.o: ../Source/paru_factorize.cpp
	$(C) -c $<

paru_fs_factorize.o: ../Source/paru_fs_factorize.cpp
	$(C) -c $<

paru_create_element.o: ../Source/paru_create_element.cpp
	$(C) -c $<

paru_assemble_row2U.o: ../Source/paru_assemble_row2U.cpp
	$(C) -c $<

paru_trsm.o: ../Source/paru_trsm.cpp
	$(C) -c $<

paru_dgemm.o: ../Source/paru_dgemm.cpp
	$(C) -c $<

paru_print.o: ../Source/paru_print.cpp
	$(C) -c $<

paru_init_rel.o: ../Source/paru_init_rel.cpp
	$(C) -c $<

paru_write.o: ../Source/paru_write.cpp
	$(C) -c $<

paru_update_rowDeg.o: ../Source/paru_update_rowDeg.cpp
	$(C) -c $<

paru_update_rel_ind.o: ../Source/paru_update_rel_ind.cpp
	$(C) -c $<


paru_cumsum.o: ../Source/paru_cumsum.cpp
	$(C) -c $<

paru_intersection.o: ../Source/paru_intersection.cpp
	$(C) -c $<

paru_bin_search.o: ../Source/paru_bin_search.cpp
	$(C) -c $<

paru_prior_assemble.o:  ../Source/paru_prior_assemble.cpp
	$(C) -c $<

paru_assemble.o: ../Source/paru_assemble.cpp
	$(C) -c $<

paru_perm.o: ../Source/paru_perm.cpp
	$(C) -c $<

paru_lsolve.o: ../Source/paru_lsolve.cpp
	$(C) -c $<

paru_usolve.o: ../Source/paru_usolve.cpp
	$(C) -c $<

paru_solve.o: ../Source/paru_solve.cpp
	$(C) -c $<

paru_residual.o: ../Source/paru_residual.cpp
	$(C) -c $<

paru_backward.o: ../Source/paru_backward.cpp
	$(C) -c $<

paru_gaxpy.o: ../Source/paru_gaxpy.cpp
	$(C) -c $<

paru_norms.o: ../Source/paru_norms.cpp
	$(C) -c $<

paru_full_summed.o: ../Source/paru_full_summed.cpp
	$(C) -c $<


#-------------------------------------------------------------------------------

# install PARU
install: $(AR_TARGET) $(INSTALL_LIB)/$(SO_TARGET)

$(INSTALL_LIB)/$(SO_TARGET): $(OBJ)
	@mkdir -p $(INSTALL_LIB)
	@mkdir -p $(INSTALL_INCLUDE)
	@mkdir -p $(INSTALL_DOC)
	$(C) $(SO_OPTS) $^ -o $@ $(LDLIBS)
	( cd $(INSTALL_LIB) ; ln -sf $(SO_TARGET) $(SO_PLAIN) )
	( cd $(INSTALL_LIB) ; ln -sf $(SO_TARGET) $(SO_MAIN) )
	$(CP) ../Include/ParU.hpp $(INSTALL_INCLUDE)
	$(CP) ../Doc/paru_user_guide.pdf $(INSTALL_DOC)
	$(CP) ../README.txt $(INSTALL_DOC)/PARU_README.txt
	chmod 755 $(INSTALL_LIB)/$(SO_TARGET)
	chmod 644 $(INSTALL_INCLUDE)/ParU.hpp
#	chmod 644 $(INSTALL_DOC)/paru_user_guide.pdf
	chmod 744 $(INSTALL_DOC)/paru_user_guide.pdf
	chmod 644 $(INSTALL_DOC)/PARU_README.txt

# uninstall PARU
uninstall:
	$(RM) $(INSTALL_LIB)/$(SO_TARGET)
	$(RM) $(INSTALL_LIB)/$(SO_PLAIN)
	$(RM) $(INSTALL_LIB)/$(SO_MAIN)
	$(RM) $(INSTALL_INCLUDE)/ParU.hpp
	$(RM) $(INSTALL_DOC)/paru_user_guide.pdf
	$(RM) $(INSTALL_DOC)/PARU_README.txt
