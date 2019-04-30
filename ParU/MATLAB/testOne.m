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

err = 10e-9;
s = 1;

%test
id = 1222;


Prob = ssget(id);
A = Prob.A;
group = index.Group {id} ;
name = index.Name {id} ;
[dp,dq,dr,ds,dcc,drr] = dmperm(A);

for s = [0 1]
    str1 = 'tar zvfxO ~/SuiteSparseCollection//MM/';
    str2 = sprintf ('%s/%s.tar.gz %s/%s.mtx | ../Demo/testazny %d %d', ...
    group, name, name, name, id, s) ;
    str = strcat (str1,str2);

    system(str);

    %%scaling now it is in the code
    %A = sparse(diag(1./max(abs(A),[],2)))*A;
    %mmwrite('../Matrix/ParUTst/tmp.mtx', A);
    %str = sprintf ('../Demo/testazny %d < ../Matrix/ParUTst/tmp.mtx', id );
    %system(str);


    % Loading the results into Matlab
    path = '../Demo/Res/';

    row_name = sprintf ('%d_row.txt', id);
    rowfullname = strcat(path, row_name);
    rowp = load (rowfullname);
    rowp = rowp+1;

    col_name = sprintf ('%d_col.txt', id);
    colfullname = strcat(path, col_name);
    colp = load (colfullname);
    colp = colp+1;

    if(s==1)
        s_name = sprintf ('%d_scale.txt', id);
        scalefullname = strcat(path, s_name);
        scale= load (scalefullname);
    end
    
    
    LU_name = sprintf ('%d_LU.txt', id);
    LUfullname = strcat(path, LU_name);
    [LU, paddingZ] = mmread (LUfullname);
    
    info_name = sprintf ('%d_info.txt',id);
    infofullname = strcat(path, info_name);
    myElaps = load (infofullname)
    
    
    
    umfStart= tic;
    %[l, u, p, q, D]=lu(A, 'vector');
    [l, u, p, q, r, Info]= umfpack (A);
    umfElaps = toc(umfStart)
    
    
    L=tril(LU,-1)+speye(size(LU));
    U=triu(LU); 
    
    if (s == 1)
        sA = sparse(diag(scale))*A;
    else
        sA = A;
    end
    
    myErr = lu_normest(sA(rowp,colp),L,U)
    %umfErr = lu_normest(D(:,p)\A(:,q),l,u)
    umfErr = lu_normest(p*(r\A)*q,l,u)
    
    umfpnnz = nnz(l)+nnz(u) - size(A,1)
    mynnz = nnz(LU) + nnz(paddingZ)
    
    %setting up KL
    opts.tol = 0; opts.btf = 0; opts.ordering = 2;
    
    
    B = sA(rowp,colp); B = B + 5*speye(size(B));
    [myx, myinfo, c]  = klu(B, opts);
    %myflop = luflop(L,U);
    myflop = myinfo.flops
    
    
    B = p*(r\A)*q; B = B + 5*speye(size(B));
    [umpfx, umpfinfo, c]  = klu(B, opts);
    %umfflop = luflop(l,u);
    umfflop = umpfinfo.flops
    
    myflop/umfflop
    
    
    if(myErr <= 100*umfErr || myErr < err)
        fprintf('Pass\n')
    else
        fprintf('Fail\n')
    end
    
    % cleaning the files because of the memory problem
    str = ['rm  ' path LU_name];    system(str);
    str = ['rm  ' path col_name];    system(str);
    str = ['rm  ' path row_name];    system(str);
    str = ['rm  ' path info_name];    system(str);
    if (s == 1)
        str = ['rm  ' path s_name];    system(str);
    end
end
