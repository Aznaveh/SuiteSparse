all:
	g++  -O3 -fexceptions -fPIC -fopenmp   -I../../include  ./paru_sym_analysis.cpp -o testazny -L/users/aznaveh/Sparse/SuiteSparse/lib -L../../lib -lspqr -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lmetis  -lm -lrt -Wl,-rpath=/users/aznaveh/Sparse/SuiteSparse/lib -llapack -lopenblas 
run:
#	./testazny < ../Matrix/c2.mtx
#	./testazny < ../Matrix/Franz6_id1959_aug.mtx
#	./testazny < ../Matrix/GD01_b.mtx
#	./testazny < ../Matrix/GD98_a.mtx
#	./testazny < ../Matrix/Groebner_id2003_aug.mtx
#	./testazny < ../Matrix/LFAT5.mtx
#	./testazny < ../Matrix/Ragusa16.mtx
	./testazny < ../Matrix/Tina_AskCal.mtx
#	./testazny < ../Matrix/a0.mtx
#	./testazny < ../Matrix/a04.mtx
#	./testazny < ../Matrix/a1.mtx
#	./testazny < ../Matrix/a4.mtx
#	./testazny < ../Matrix/a2.mtx
#	./testazny < ../Matrix/arrow.mtx
#	./testazny < ../Matrix/ash219.mtx
	./testazny < ../Matrix/b1_ss.mtx
#	./testazny < ../Matrix/test.mtx

