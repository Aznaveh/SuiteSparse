
clear

index = ssget ;

f = find (index.nrows == index.ncols & ...
    index.sprank == index.ncols & ...
    ~index.posdef & ...
    index.isReal & ~index.isGraph) ;

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
        fprintf ('%20s/%-30s n: %10d nnz: %10d\n', ...
            index.Group {id}, index.Name {id}, ...
            index.nrows (id), thisnz) ;
    end
    nzlast = thisnz ;
end

fnew = fnew' ;
nmat = length (fnew) ;

% tar zvfxO ~/ssget/MM/HB/west0067.tar.gz west0067/west0067.mtx | \
%     ./build/ttest

ff = fopen ('myscript', 'w') ;
for k = 1:nmat
    id = fnew (k) ;
    group = index.Group {id} ;
    name = index.Name {id} ;
    fprintf (ff,...
    '\ntar zvfxO ~/SuiteSparseCollection//MM/%s/%s.tar.gz %s/%s.mtx | ./testazny', ...
       group, name, name, name) ;
end
fclose (ff) ;
