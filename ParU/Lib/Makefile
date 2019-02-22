#===============================================================================
# ParU/Lib/Makefile
#===============================================================================

LIBRARY = libparu
VERSION = 0.0.1
SO_VERSION = 0

default: library

ccode: all

include ../../SuiteSparse_config/SuiteSparse_config.mk


# ParU depends on SPQR, CHOLMOD, AMD, COLAMD, LAPACK, the BLAS and SuiteSparse_config
LDLIBS += -lspqr -lamd -lcolamd -lcholmod -lsuitesparseconfig $(LAPACK) $(BLAS)

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

INC = ../Include/Parallel_LU.hpp

OBJ = \
    paru_mem.o \
    paru_sym_analyse.o \
	paru_init_rowFronts.o\
	paru_tuples.o\
	paru_assemble.o\
	paru_factorize.o\
	paru_create_element.o\
	paru_numeric_assemble.o\
	paru_trsm.o\
	paru_dgemm.o\
	paru_print.o\
	paru_fourPass.o\
	paru_init_rel.o\
	paru_update_rel_ind.o

$(OBJ): $(INC)

I = -I../../SPQR/Include -I../../CHOLMOD/Include -I../../SuiteSparse_config -I../Include

C = $(CXX) $(CF) $(SPQR_CONFIG) $(CHOLMOD_CONFIG) $(CONFIG_PARTITION) $(I)
#-------------------------------------------------------------------------------

static: $(AR_TARGET)

$(AR_TARGET): $(OBJ)
	$(ARCHIVE) $@ $^
	- $(RANLIB) $@

#-------------------------------------------------------------------------------

paru_mem.o: ../Source/paru_mem.cpp
	$(C) -c $<

paru_sym_analyse.o: ../Source/paru_sym_analyse.cpp
	$(C) -c $<

paru_init_rowFronts.o: ../Source/paru_init_rowFronts.cpp
	$(C) -c $<

paru_tuples.o: ../Source/paru_tuples.cpp
	$(C) -c $<

paru_assemble.o: ../Source/paru_assemble.cpp
	$(C) -c $<

paru_factorize.o: ../Source/paru_factorize.cpp
	$(C) -c $<

paru_create_element.o: ../Source/paru_create_element.cpp
	$(C) -c $<

paru_numeric_assemble.o: ../Source/paru_numeric_assemble.cpp
	$(C) -c $<

paru_trsm.o: ../Source/paru_trsm.cpp
	$(C) -c $<

paru_dgemm.o: ../Source/paru_dgemm.cpp
	$(C) -c $<

paru_fourPass.o: ../Source/paru_fourPass.cpp
	$(C) -c $<

paru_print.o: ../Source/paru_print.cpp
	$(C) -c $<

paru_init_rel.o: ../Source/paru_init_rel.cpp
	$(C) -c $<

paru_update_rel_ind.o: ../Source/paru_update_rel_ind.cpp
	$(C) -c $<




#-------------------------------------------------------------------------------

# install PARU
install: $(AR_TARGET) $(INSTALL_LIB)/$(SO_TARGET)

$(INSTALL_LIB)/$(SO_TARGET): $(OBJ)
	@mkdir -p $(INSTALL_LIB)
	@mkdir -p $(INSTALL_INCLUDE)
	@mkdir -p $(INSTALL_DOC)
	$(CXX) $(SO_OPTS) $^ -o $@ $(LDLIBS)
	( cd $(INSTALL_LIB) ; ln -sf $(SO_TARGET) $(SO_PLAIN) )
	( cd $(INSTALL_LIB) ; ln -sf $(SO_TARGET) $(SO_MAIN) )
	$(CP) ../Include/Parallel_LU.hpp $(INSTALL_INCLUDE)
	$(CP) ../Doc/paru_user_guide.pdf $(INSTALL_DOC)
	$(CP) ../README.txt $(INSTALL_DOC)/PARU_README.txt
	chmod 755 $(INSTALL_LIB)/$(SO_TARGET)
	chmod 644 $(INSTALL_INCLUDE)/Parallel_LU.hpp
#	chmod 644 $(INSTALL_DOC)/paru_user_guide.pdf
	chmod 744 $(INSTALL_DOC)/paru_user_guide.pdf
	chmod 644 $(INSTALL_DOC)/PARU_README.txt

# uninstall PARU
uninstall:
	$(RM) $(INSTALL_LIB)/$(SO_TARGET)
	$(RM) $(INSTALL_LIB)/$(SO_PLAIN)
	$(RM) $(INSTALL_LIB)/$(SO_MAIN)
	$(RM) $(INSTALL_INCLUDE)/Parallel_LU.hpp
	$(RM) $(INSTALL_DOC)/paru_user_guide.pdf
	$(RM) $(INSTALL_DOC)/PARU_README.txt

