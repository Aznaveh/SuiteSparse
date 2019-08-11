% A must be there 
id = 0
s = 0
[m n] = size (A);
err = 10e-9;
mmwrite('../Matrix/ParUTst/tmp.mtx', A);
%str = sprintf ('../Demo/testazny %d %d < ../Matrix/ParUTst/tmp.mtx', id, s);
str = sprintf ('../Demo/umfout %d %d < ../Matrix/ParUTst/tmp.mtx', id, s);
system(str);


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
    Rvec = load (scalefullname);
    R = spdiags (Rvec, 0, m, m);
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
    %sA = sparse(diag(scale))*A;
    sA = R \A;
else
    sA = A;
end

myErr = lu_normest(sA(rowp,colp),L,U)
%umfErr = lu_normest(D(:,p)\A(:,q),l,u)
umfErr = lu_normest(p*(r\A)*q,l,u)

umfpnnz = nnz(l)+nnz(u) - m;
mynnz = nnz(LU); %+ nnz(paddingZ)

%setting up KL
%opts.tol = 0; opts.btf = 0; opts.ordering = 2;
%B = sA(rowp,colp); B = B + 5*speye(size(B));
%[myx, myinfo, c]  = klu(B, opts);
%myflop = myinfo.flops
myflop = luflop(L,U);


%B = p*(r\A)*q; B = B + 5*speye(size(B));
%[umpfx, umpfinfo, c]  = klu(B, opts);
%umfflop = umpfinfo.flops
umfflop = luflop(l,u);

myflop/umfflop;


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
