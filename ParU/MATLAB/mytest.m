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
       % fprintf ('%20s/%-30s n: %10d nnz: %10d\n', ...
       %     index.Group {id}, index.Name {id}, ...
       %     index.nrows (id), thisnz) ;
    end
    nzlast = thisnz ;
end

fnew = fnew' ;
nmat = length (fnew) ;


err = 1e-5;

ff = fopen ('results.out', 'w') ;

% Headers
fprintf(ff,'id\tmyErr\tumfErr\tratio' );
fprintf(ff,'\tmyElaps\tumfElaps\tratio');
fprintf(ff,'\tmynnz\tumfnnz\tratio');
fprintf(ff,'\tmyflop\tumfflop\tratio\n');


for k = 1:nmat
    id = fnew (k) ;
    group = index.Group {id} ;
    name = index.Name {id} ;

    Prob = ssget(id);
    A = Prob.A;

    str1 = 'tar zvfxO ~/SuiteSparseCollection//MM/';
    str2 = sprintf ('%s/%s.tar.gz %s/%s.mtx | ../Demo/testazny %d', ...
    group, name, name, name, id) ;
    str = strcat (str1,str2);

    myStart = tic;
    system(str);
    myElaps = toc(myStart);


    % Loading the results into Matlab
    %path = '/export/scratch/multifrontal/aznaveh/myResults/';
    path = '../Demo/Res/';

    row_name = sprintf ('%d_row.txt',id);
    rowfullname = strcat(path, row_name);
    rowp = load (rowfullname);
    rowp = rowp+1;

    col_name = sprintf ('%d_col.txt',id);
    colfullname = strcat(path, col_name);
    colp = load (colfullname);
    colp = colp+1;

    LU_name = sprintf ('%d_LU.txt',id);
    LUfullname = strcat(path, LU_name);

    [LU, paddingZ] = mmread (LUfullname);


    umfStart= tic;
    %[l,u,p,q,D]=lu(A, 'vector');
    [l, u, p, q, r, Info]= umfpack (A);
    umfElaps = toc(umfStart);


    L=tril(LU,-1)+speye(size(LU));
    U=triu(LU); 

    myErr = lu_normest(A(rowp,colp),L,U);
    %umfErr = lu_normest(D(:,p)\A(:,q),l,u);
    umfErr = lu_normest(p*(r\A)*q,l,u);

    nnzumfp = nnz(l)+nnz(u) - size(A,1);
    mynnz = nnz(LU) + nnz(paddingZ);
    myflop = luflop(L,U);
    umfflop = luflop(l,u);

    fprintf(ff,'%d\t%g\t%g\t%g', id, myErr, umfErr, myErr/umfErr);
    fprintf(ff,'\t%g\t%g\t%g', myElaps, umfElaps, myElaps/umfElaps);
    fprintf(ff,'\t%g\t%g\t%g', mynnz , nnzumfp, mynnz/nnzumfp );
    fprintf(ff,'\t%g\t%g\t%g', myflop, umfflop, myflop/umfflop);


    if(myErr <= 100*umfErr || myErr < err)
        fprintf(ff,'\tPass\n');
    else
        fprintf(ff,'\tFail\n');
    end

    % cleaning the files because of the memory problem
    str = ['rm  ' path LU_name];    system(str);
    str = ['rm  ' path col_name];    system(str);
    str = ['rm  ' path row_name];    system(str);
end
fclose (ff) ;
