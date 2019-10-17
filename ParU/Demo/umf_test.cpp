/** =========================================================================  /
 * =======================  paru_test========================================  /
 * ==========================================================================  /
 * @brief    test to see how to call umfpack symbolic analysis
 * 
 * @author Aznaveh
 * */
#include "Parallel_LU.hpp"
#include <omp.h>

int main (int argc, char **argv)
{
    DEBUGLEVEL(0); 

    cholmod_common Common, *cc;
    cholmod_sparse *A;
    int mtype;
    paru_symbolic *LUsym;


    //~~~~~~~~~Reading the input matrix and test if the format is OK~~~~~~~~~~~~
    // start CHOLMOD
    cc = &Common;
    cholmod_l_start (cc);

    // A = mread (stdin) ; read in the sparse matrix A
    A = (cholmod_sparse *) cholmod_l_read_matrix (stdin, 1, &mtype, cc);
    if (A == NULL){
        printf ("Paru: input matrix is invalid\n");
        exit (1);
    }

    if (mtype != CHOLMOD_SPARSE){
        printf ("Paru: input matrix must be sparse\n");
        exit (1);
    }

    if (A->xtype != CHOLMOD_REAL){
        printf ("Paru: input matrix must be real\n");
        exit (1);
    }

    //~~~~~~~~~~~~~~~~~~~Starting computation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double my_start_time = omp_get_wtime();

    LUsym = paru_analyze (A, cc);
    if (LUsym == NULL) {
        exit(0);
    }


    int scale = 0;
    if (argc == 3){
        scale = atoi(argv[2]);
        if (scale) 
            PRLEVEL (1, ("The input matrix will be scaled\n"));
    }
        
    paru_matrix *paruMatInfo = paru_init_rowFronts (A, scale, LUsym, cc);
    if (paruMatInfo == NULL) {
        exit(0);
    }


    //Int m = paruMatInfo-> m;
    //Int n = paruMatInfo-> n;
    Int nf = paruMatInfo->LUsym->nf;


    for (Int i = 0; i < nf; i++) {
        if (paru_front (paruMatInfo, i, cc)){
            printf ("some problem\n");
            break;
        }
    }
    double my_time = omp_get_wtime() - my_start_time;
    paruMatInfo->my_time = my_time;
 
    //~~~~~~~~~~~~~~~~~~~End computation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    //~~~~~~~~~~~~~~~~~~~Calling umfpack~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    double status,   // Info [UMFPACK_STATUS] 
    Info[UMFPACK_INFO],// Contains statistics about the symbolic analysis

    Control [UMFPACK_CONTROL]; // it is set in umfpack_dl_defaults and
    // is used in umfpack_dl_symbolic; if
    // passed NULL it will use the defaults
    umfpack_dl_defaults (Control) ;
    Control [UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;

    Int *Ap = (Int*) A->p;
    Int *Ai = (Int*) A->i;
    double *Ax = (double*) A->x;
    Int m = A->nrow;
    Int n = A->ncol;
    void *Symbolic, *Numeric;  // Output argument in umf_dl_symbolc;

    double umf_start_time = omp_get_wtime();
    status = umfpack_dl_symbolic (n, n, Ap, Ai, Ax, &Symbolic,
            Control, Info) ;
    if (status < 0)
    {
        umfpack_dl_report_info (Control, Info) ;
        umfpack_dl_report_status (Control, status) ;
        printf ("umfpack_dl_symbolic failed\n") ;
    }
    status = umfpack_dl_numeric (Ap, Ai, Ax, Symbolic, &Numeric,
	Control, Info) ;
    if (status < 0)
    {
	umfpack_dl_report_info (Control, Info) ;
	umfpack_dl_report_status (Control, status) ;
	printf ("umfpack_dl_numeric failed\n") ;
    }

    double umf_time = omp_get_wtime() - umf_start_time;
    umfpack_dl_free_symbolic (&Symbolic) ;
    umfpack_dl_free_numeric (&Numeric) ;

    paruMatInfo->umf_time = umf_time;



    // Writing results to a file
    paru_write(paruMatInfo, scale,  argv[1], cc);

    //~~~~~~~~~~~~~~~~~~~Free Everything~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cholmod_l_free_sparse (&A, cc);

    paru_freemat (&paruMatInfo, cc);
    paru_freesym (&LUsym,cc);

    cholmod_l_finish (cc);

}
