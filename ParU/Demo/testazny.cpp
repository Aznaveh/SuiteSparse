/** =========================================================================  /
 * =======================  paru_test========================================  /
 * ==========================================================================  /
 * @brief   doing the LU factorization of the matrix one by one
 *          argv[1] = name of the input matrix that is used for output
 *                  default is zero
 *          argv[2] = 0 if it doesnt needs to be scaled
 *                  default is no scaling
 *
 * @author Aznaveh
 * */
#include "Parallel_LU.hpp"
#include <omp.h>

#define PRINTCBsTUPLES 0
int main (int argc, char **argv)
{
    DEBUGLEVEL(0); 
    PRLEVEL (0, ("%%%%%%%%%%%%%%%%%% START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n"));
    PRLEVEL (0, ("clear all\n" ));
    cholmod_common Common, *cc;
    cholmod_sparse *A;
    int mtype;
    paru_symbolic *LUsym;

    char name[25];
    //fgets(name, 25, stdin);
    //fscanf (stdin, "%s", &name);
    //printf("The name is %s\n", name);
    //fseek(stdin, 0, SEEK_SET);

    


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


    LUsym = paru_sym_analyse (A, cc);
    if (LUsym == NULL) {
        exit(0);
    }


    //~~~~~~~~~~~~~~~~~~~Starting computation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double start_time = omp_get_wtime();

    int scale = 0;
    if (argc == 3){
        scale = atoi(argv[2]);
        if (scale) printf ("The input matrix will be scaled\n");
    }
        
    paru_matrix *paruMatInfo = paru_init_rowFronts (A, scale, LUsym, cc);
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
        //paru_assemble (paruMatInfo, i, cc);
        paru_front (paruMatInfo, i, cc);
    }
    double time = omp_get_wtime() - start_time;
    paruMatInfo->time = time;
    //~~~~~~~~~~~~~~~~~~~End computation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    paru_write(paruMatInfo, scale,  argv[1], cc);

    cholmod_l_free_sparse (&A, cc);

    paru_freemat (&paruMatInfo, cc);
    paru_freesym (&LUsym,cc);

    cholmod_l_finish (cc);

    //Matlab
    //
    Int p = 0;
    PRLEVEL (p, ("err = 1e-9; [m n] =size(S);\n"));
    PRLEVEL (p, ("format long g\n"));
    PRLEVEL (p, ("oldR=[]; c=[];\n"));
    PRLEVEL (p, ("%%Finalizing the permutation\n"));
    PRLEVEL (p, ("for f=1:%ld\n",nf));
    PRLEVEL (p, ("\tnpivots(f) = length(cols{f});\n"));
    PRLEVEL (p, (
                "\toldR=[oldR, rows{f}(1:npivots(f))]; c=[c, cols{f}];\nend\n"));
    PRLEVEL (p, ("oldR = [oldR setdiff(1:m,oldR)];\n"));
    PRLEVEL (p, ("newR(oldR)=1:length(oldR);\n"));

    PRLEVEL (p, ("for f=1:%ld\n",nf));
    PRLEVEL (p, ("\tLU(newR(rows{f}),cols{f})=Luf{f};\n"));
    PRLEVEL (p, ("\tLU(newR(Urows{f}),Ucols{f})=Us{f};\nend\n"));

    PRLEVEL (p, ("L=tril(LU,-1)+eye(size(LU));\n"));
    PRLEVEL (p, ("U=triu(LU); U=U(1:n,:);\n"));
    PRLEVEL (p, ("U=sparse(U); L=sparse(L);\n"));
    PRLEVEL (p, ("spparms('spumoni',3);\n" ));
    PRLEVEL (p, ("[l,u,p]=lu(S, 'vector');\n" ));
    //PRLEVEL (p, ("matlabErr = norm(S(p,:)-l*u,'fro')\n" ));
    PRLEVEL (p, ("matlabErr = lu_normest(S(p,:),l,u)\n" ));
    PRLEVEL (p, ("nnzMat= nnz(l)+nnz(u) \n" ));
    PRLEVEL (p, ("fprintf('Paru\\n');\n"));
    //PRLEVEL (p, ("myErr = norm(S(oldR,c)-L*U,'fro')\n" ));
    PRLEVEL (p, ("myErr = lu_normest(S(oldR,c),L,U)\n" ));
    PRLEVEL (p, ("mynnz= nnz(L)+nnz(U) \n" ));
    //PRLEVEL (p, ("if( (norm(S(oldR,c)-L*U)) < err )\n" ));
    PRLEVEL (p, ("if(myErr <= 100*matlabErr || myErr<err)\n" ));
    PRLEVEL (p, ("\tfprintf('Pass\\n')\nelse\n\tfprintf('Fail\\n')\nend\n" ));
    PRLEVEL (p, ("pause\n" ));

    PRLEVEL (p, ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"));
    //   printf 
    //   ("malloc_count %ld inuse %ld\n", cc->malloc_count, cc->memory_inuse);
}
