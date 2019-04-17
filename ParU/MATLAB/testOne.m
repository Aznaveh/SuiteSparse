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


id = 1576; 
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
InMatrix = load (LUfullname);


I = InMatrix(:,1)+1;
J = InMatrix(:,2)+1;
X = InMatrix(:,3);
LU = sparse(I,J,X);


%umfStart= tic;
%[l,u,p]=lu(A, 'vector');
%umfElaps = toc(umfStart);


L=tril(LU,-1)+speye(size(LU));
U=triu(LU); 
myErr = lu_normest(A(rowp,colp),L,U);
%matlabErr = lu_normest(A(p,:),l,u);

if(myErr <= 100*matlabErr || myErr < err)
    fprintf('Pass\n')
else
    fprintf('Fail\n')
end
