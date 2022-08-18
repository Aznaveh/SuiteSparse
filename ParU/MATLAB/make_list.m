ff = fopen ('list.sh', 'w') ;
coll_path = '/home/grads/a/aznaveh/SuiteSparse/ssget/MM';
coll_path = '/raid/archive/davis/SuiteSparseCollection/MM';
fprintf(ff,'export COL_PATH=$COL_PATH"%s"\n',coll_path );
fprintf(ff,'echo "Starting a suite of tests (the list is generated by MATLAB)"\n');

index = ssget ;

f = find (index.nrows == index.ncols & ...
index.sprank == index.ncols & ...
index.sprank == index.nrows & ...
index.isReal & ~index.isGraph & ...
index.numerical_symmetry < .15 &...
~index.posdef & ...
index.nnz >=1e6 & ...
index.pattern_symmetry <= .65 ) ;

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
    fprintf(ff,'tar zvfxO $COL_PATH/%s/%s.tar.gz',group,name);
    fprintf(ff,' | ./x_paru_demo\n');
end
