/** =========================================================================  /
 * =======================  matlab_finalize =================================  /
 * ==========================================================================  /
 * @brief    To use the matlab out put
 * 
 * @author Aznaveh
 * */

#include "Parallel_LU.hpp"
void matlab_finalize (Int nf)
{

    Int p=0;
    PRLEVEL (p, ("err = 1e-9; [m n] =size(S);\n"));
    PRLEVEL (p, ("format long g\n"));
    PRLEVEL (p, ("oldR=[]; c=[];\n"));
    PRLEVEL (p, ("%%Finalizing the permutation\n"));
    PRLEVEL (p, ("for f=1:%ld\n",nf));
    PRLEVEL (p, ("\tnpivots(f) = length(cols{f});\n"));
    PRLEVEL (p, ("\toldR=[oldR, rows{f}(1:npivots(f))]; c=[c, cols{f}];\nend\n"));
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
//    printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
}
