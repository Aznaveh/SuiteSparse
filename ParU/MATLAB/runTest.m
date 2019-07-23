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
    end
    nzlast = thisnz ;
end

fnew = fnew' ;
nmat = length (fnew) ;


err = 1e-5;
%Don't scale the matrix if s ==0 scale otherwise
s = 0;

ff = fopen ('myRes.m', 'w') ;

% Headers
fprintf(ff,'id nnzA myErr umfErr logratio' );
fprintf(ff,' myElaps umfElaps ratio');
fprintf(ff,' mynnz umfnnz ratio');
fprintf(ff,' myflop umfflop ratio\n');


for k = 1:nmat
%for k = 1:10
    id = fnew (k) ;
    group = index.Group {id} ;
    name = index.Name {id} ;

    Prob = ssget(id);
    A = Prob.A;
    [dp,dq,dr,ds,dcc,drr] = dmperm(A);

    [m n] = size (A);
    if (size(dr) ~= 2 )
       continue 
       % if (norm(diff(dr)-diff(ds)) ~= 0 )
       %     sprintf('Unexpected')
       %     continue;
       % end
       % B = A(dp,dq);
       % [M,I] = max(diff(dr));
       % A = B(dr(I):dr(I+1)-1, dr(I):dr(I+1)-1 );

       % mmwrite('../Matrix/ParUTst/tmp.mtx', A);
       % str = sprintf ('../Demo/umfout %d < ../Matrix/ParUTst/tmp.mtx', id );
       % system(str);
    else 
        str1 = 'tar zvfxO ~/SuiteSparseCollection//MM/';
        str2 = sprintf ('%s/%s.tar.gz %s/%s.mtx | ../Demo/umfout %d %d', ...
        group, name, name, name, id, s) ;
        dostr = strcat (str1,str2);
        system(dostr);

        %scaling
        %A = sparse(diag(1./max(abs(A),[],2)))*A;
        %mmwrite('../Matrix/ParUTst/tmp.mtx', A);
        %str = sprintf ('../Demo/testazny %d < ../Matrix/ParUTst/tmp.mtx', id );
        %system(str);

    end


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

    if(s == 1)
        s_name = sprintf ('%d_scale.txt', id);
        scalefullname = strcat(path, s_name);
        Rvec = load (scalefullname);
        R = spdiags (Rvec, 0, m, m);
    end


    LU_name = sprintf ('%d_LU.txt',id);
    LUfullname = strcat(path, LU_name);

    info_name = sprintf ('%d_info.txt',id);
    infofullname = strcat(path, info_name);
    myElaps = load (infofullname);


    [LU, paddingZ] = mmread (LUfullname);


    umfStart= tic;
    %[l,u,p,q,D]=lu(A, 'vector');
    [l, u, p, q, r, Info]= umfpack (A);
    umfElaps = toc(umfStart);


    L=tril(LU,-1)+speye(size(LU));
    U=triu(LU); 

    if (s == 1)
        sA = R \A;
    else
        sA = A;
    end
    myErr = lu_normest(sA(rowp,colp),L,U);
    %umfErr = lu_normest(D(:,p)\A(:,q),l,u);
    umfErr = lu_normest(p*(r\A)*q,l,u);

    umfpnnz= nnz(l)+nnz(u) - m;
    mynnz = nnz(LU); %+ nnz(paddingZ);



    myflop = luflop(L,U);

    umfflop = luflop(l,u);

    fprintf(ff,'%d %d %g %g %g', id, nnz(A), myErr, umfErr, ...
        log10(myErr/umfErr));
    fprintf(ff,' %g %g %g', myElaps, umfElaps, myElaps/umfElaps);
    fprintf(ff,' %g %g %g', mynnz , umfpnnz, mynnz/umfpnnz );
    fprintf(ff,' %g %g %g\n', myflop, umfflop, myflop/umfflop);


    % cleaning the files because of the memory problem
    str = ['rm  ' path LU_name];    system(str);
    str = ['rm  ' path col_name];    system(str);
    str = ['rm  ' path row_name];    system(str);
    str = ['rm  ' path info_name];    system(str);
    if (s == 1)
        str = ['rm  ' path s_name];    system(str);
    end
end
fclose (ff) ;
