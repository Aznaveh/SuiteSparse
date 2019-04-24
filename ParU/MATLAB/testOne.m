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

% Bad Error
id = 956; 
id = 2269; 
id = 350;
id = 212; 

% Good one vs UMFPACK
id = 113; 

%test
id = 1294; 

Prob = ssget(id);
A = Prob.A;
group = index.Group {id} ;
name = index.Name {id} ;


str1 = 'tar zvfxO ~/SuiteSparseCollection//MM/';
str2 = sprintf ('%s/%s.tar.gz %s/%s.mtx | ../Demo/testazny %d', ...
group, name, name, name, id) ;
str = strcat (str1,str2);

myStart = tic;
system(str);
myElaps = toc(myStart)


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

LU_name = sprintf ('%d_LU.txt', id);
LUfullname = strcat(path, LU_name);
[LU, paddingZ] = mmread (LUfullname);

%	InMatrix = load (LUfullname);
%	I = InMatrix(:,1)+1;
%	J = InMatrix(:,2)+1;
%	X = InMatrix(:,3);
%	LU = sparse(I,J,X);



umfStart= tic;
%[l, u, p, q, D]=lu(A, 'vector');
[l, u, p, q, r, Info]= umfpack (A);
umfElaps = toc(umfStart)


L=tril(LU,-1)+speye(size(LU));
U=triu(LU); 
myErr = lu_normest(A(rowp,colp),L,U)
%umfErr = lu_normest(D(:,p)\A(:,q),l,u)
umfErr = lu_normest(p*(r\A)*q,l,u)

nnzumfp = nnz(l)+nnz(u) - size(A,1)
mynnz = nnz(LU) + nnz(paddingZ)
myflop = luflop(L,U)
umfflop = luflop(l,u)



if(myErr <= 100*umfErr || myErr < err)
    fprintf('Pass\n')
else
    fprintf('Fail\n')
end
