clear

index = ssget ;

f = find (index.nrows == index.ncols & ...
index.sprank == index.ncols & ...
~index.posdef & ...
index.isReal & ~index.isGraph & ...
index.pattern_symmetry <= .5 ) ;

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


err = 1e-5;
%Don't scale the matrix if s ==0 scale otherwise
s = 0;

ff = fopen ('list.txt', 'w') ;

% Headers

fprintf(ff,'%% id nnzA \n' );


for k = 1:nmat
    id = fnew (k) 
    % some problem in these matrice
    group = index.Group {id} ;
    name = index.Name {id} ;

    Prob = ssget(id);
    A = Prob.A;


    [dp,dq,dr,ds,dcc,drr] = dmperm(A);

    [m n] = size (A);
    if (size(dr) ~= 2 )
        %continue 
        if (norm(diff(dr)-diff(ds)) ~= 0 )
            sprintf('Unexpected')
            continue;
        end
        B = A(dp,dq);
        [M,I] = max(diff(dr));
        A = B(dr(I):dr(I+1)-1, dr(I):dr(I+1)-1 );


        [m n] = size (A);
        if ( m== 1)
            sprintf('not worth trying');
            continue;
        end
    end

    if (nnz(A) < 1000)
            continue;
    end;

    fprintf(ff,'%d %d \n', id, nnz(A)); 
end

fclose (ff) ;
