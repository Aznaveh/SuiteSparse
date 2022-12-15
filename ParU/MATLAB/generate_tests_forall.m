paru = fopen ('paru.sh', 'w') ;
mumps = fopen ('mumps.sh', 'w') ;
superlu = fopen ('superlu.sh', 'w') ;

coll_path = '/home/grads/a/aznaveh/SuiteSparse/ssget/MM';
coll_path = '/raid/archive/davis/SuiteSparseCollection/MM';

fprintf(paru,'export COL_PATH=$COL_PATH"%s"\n',coll_path );
fprintf(paru,'echo "Starting a suite of tests (the list is generated by MATLAB)"\n');

fprintf(mumps,'export COL_PATH=$COL_PATH"%s"\n',coll_path );
fprintf(mumps,'echo "Starting a suite of tests (the list is generated by MATLAB)"\n');

coll_path_superlu = '/raid/archive/davis/SuiteSparseCollection/RB';
fprintf(superlu,'export COL_PATH=$COL_PATH"%s"\n',coll_path_superlu );
fprintf(superlu,'echo "Starting a suite of tests (the list is generated by MATLAB)"\n');


index = ssget ;

f = find (index.nrows == index.ncols & ...
index.sprank == index.ncols & ...
index.sprank == index.nrows & ...
index.isReal & ~index.isGraph & ...
index.numerical_symmetry <= .15 &...
~index.posdef & ...
index.nnz >1e6 & ...
index.pattern_symmetry <= .7 ) ;

[ignore, i] = sort (index.nnz (f) + index.nzero (f)) ;

f = f (i) ;

nzlast = -1 ;
fnew = [ ] ;

nmat = length (f) ;
for k = 1:nmat
    id = f (k) ;

    thisnz = index.nnz (id) + index.nzero (id) ;
    if (thisnz ~= nzlast)
        fnew = [fnew id] ;
    end
    nzlast = thisnz ;
end

fnew = fnew' ;

nmat = length (fnew) ;


for k = 1:nmat
    id = fnew (k); 
    group = index.Group {id} ;
    name = index.Name {id} ;

    fprintf(paru,'tar zvfxO $COL_PATH/%s/%s.tar.gz',group,name);
    fprintf(paru,' | ./x_paru_demo\n');

    fprintf(mumps,'tar zvfxO $COL_PATH/%s/%s.tar.gz',group,name);
    fprintf(mumps,' | ./mumps_paru | grep ''Elapsed time in'' | tail -n 3 | awk ''NF>1{print $NF}'' \n');

    fprintf(superlu,'tar zvfxO $COL_PATH/%s/%s.tar.gz',group,name);
    fprintf(superlu,'|  ./pdlinsolx -p 24 | grep ''time'' | awk ''NF>1{print $NF}'' \n');

end