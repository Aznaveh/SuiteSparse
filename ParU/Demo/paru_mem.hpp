/* Wrappers for managing memory */
//#include "Parallel_LU.hpp"
#include"spqr.hpp"
void *paralloc(int n, int size, cholmod_common* cc);
void paru_freesym(paru_symbolic** LUsym_handle,cholmod_common *cc);

