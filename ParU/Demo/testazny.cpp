#include "Parallel_LU.hpp"
#define PRINTCBsTUPLES 0

int main (int argc, char **argv)
{
    DEBUGLEVEL(0); 
    printf("%%%%%%%%%%%%%%%%%% START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
    PRLEVEL (0, ("clear all\n" ));
    cholmod_common Common, *cc;
    cholmod_sparse *A;
    int mtype;
    paru_symbolic *LUsym;

    // start CHOLMOD
    cc = &Common;
    cholmod_l_start (cc);

    // A = mread (stdin) ; read in the sparse matrix A
    A = (cholmod_sparse *) cholmod_l_read_matrix (stdin, 1, &mtype, cc);
    if (mtype != CHOLMOD_SPARSE)
    {
        printf ("input matrix must be sparse\n");
        exit (1);
    }

    LUsym = paru_sym_analyse (A, cc);
    if (LUsym == NULL) {
        exit(0);
    }


    paru_matrix *paruMatInfo = paru_init_rowFronts (A, LUsym, cc);
    if (paruMatInfo == NULL) {
        exit(0);
    }

    Int m,n,nf;
    m = paruMatInfo-> m;
    n = paruMatInfo-> n;
    nf = paruMatInfo->LUsym->nf;

    if (PRINTCBsTUPLES){
        for (int i = 0; i < m+nf+1; ++i) 
            paru_print_element (paruMatInfo, i);

        tupleList *RowList = paruMatInfo -> RowList;
        tupleList *ColList = paruMatInfo -> ColList;  

        printf ("RowList:\n");
        for (Int i = 0; i < m; ++i){ 
            printf("row %ld :",i);
            paru_print_tupleList (RowList , i);
        }

        PRLEVEL (1, ("ColList =%p\n", ColList));
        printf ("ColList:\n");
        for (Int i = 0; i < n; ++i) {
            printf("col %ld :",i);
            paru_print_tupleList (ColList , i);
        }

    }

    for (Int i = 0; i < nf; i++) {
        paru_assemble (paruMatInfo, i, cc);
    }



    cholmod_l_free_sparse (&A, cc);
    paru_freemat (&paruMatInfo, cc);

    paru_freesym (&LUsym,cc);


    //  spqr_symbolic *QRsym = 
    //      spqr_analyze (A, SPQR_ORDERING_CHOLMOD, FALSE,FALSE , FALSE, cc);
    //           spqr_freesym (&QRsym, cc);
    //  PRLEVEL (1, ("malloc_count %ld inuse %ld\n", 
    //cc->malloc_count, cc->memory_inuse));

    cholmod_l_finish (cc);

    //Matlab
    //
    Int p=0;
    PRLEVEL (p, ("ncols=0;\n invp(p)=1:%ld;\n",n));
    PRLEVEL (p, ("for f=1:%ld\n",nf));
    PRLEVEL (p, ("oldrows = rows{f};\n "));
    PRLEVEL (p, ("newrows= invp(oldrows);\n"));
    PRLEVEL (p, ("newcols = ncols+1: ncols+npivots(f); \n"));
    PRLEVEL (p, ("\tLU(newrows,newcols)=Luf{f};\n"));
    PRLEVEL (p, ("\tLU(newrows,ncols+1:ncols+npivots(f))=U{f};\nend\n"));
    PRLEVEL (p, (" ncols = ncols+npivots(f); \n"));
    PRLEVEL (p, ("if( norm(LU-lu(S(invp,:))) < eps )\n" ));
    PRLEVEL (p, ("\tfprintf('Pass\\n')\nelse\nfprintf('Fail\\n')\nend\n" ));
    
//    PRLEVEL (p, ("U =triu(A);"));
//    PRLEVEL (p, ("d =triu(tril(A));"));
//    PRLEVEL (p, ("L =tril(A)-d+eye(%ld,%ld);\n",m,n));
    printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*\n");
    //   printf 
    //   ("malloc_count %ld inuse %ld\n", cc->malloc_count, cc->memory_inuse);
}
