%function A = testWithID(id)

err = 10e-9;
s = 0;

index = ssget ;

Prob = ssget(id);
A = Prob.A;
[m n] = size (A);
group = index.Group {id} ;
name = index.Name {id} ;
[dp,dq,dr,ds,dcc,drr] = dmperm(A);

[m n] = size (A);
if (size(dr) ~= 2 )
    if (norm(diff(dr)-diff(ds)) ~= 0 )
        sprintf('Unexpected')
    end
    B = A(dp,dq);
    [M,I] = max(diff(dr));
    A = B(dr(I):dr(I+1)-1, dr(I):dr(I+1)-1 );

    [m n] = size (A);
    if ( m== 1)
        sprintf('not worth trying');
    end

end

%max scaling
A = spdiags (1./max (A,[], 2), 0, size(A,1), size(A,2)) * A ;


mmwrite('../Matrix/ParUTst/tmp.mtx', A);
intel = sprintf('. /home/grads/a/aznaveh/intel/bin/compilervars.sh intel64;');
intel = sprintf('. /opt/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64;');
str = sprintf ('../Demo/umfout %d %d< ../Matrix/ParUTst/tmp.mtx', id, s);
%str = strcat(intel, str);

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
        Rvec = Rvec (rowp);
        R = spdiags (Rvec, 0, m, m);
    else 
        r=1;
end 


LU_name = sprintf ('%d_LU.txt', id);
LUfullname = strcat(path, LU_name);
[LU, paddingZ] = mmread (LUfullname);

info_name = sprintf ('%d_info.txt',id);
infofullname = strcat(path, info_name);
t_Info = load (infofullname);
myElaps = t_Info(1)
fromCode_umf_Elaps = t_Info(2);

flp_cnt_dgemm = t_Info(3);
flp_cnt_trsm = t_Info(4);
flp_cnt_dger = t_Info(5);
hardware_flp_cnt = flp_cnt_dgemm + flp_cnt_trsm + flp_cnt_dger;




umfStart= tic;
%[l, u, p, q, D]=lu(A, 'vector');
%[l, u, p, q, r, Info]= umfpack (A);
[l, u, p, q ]= umfpack (A);
%umfElaps = toc(umfStart)
umfElaps = fromCode_umf_Elaps


L=tril(LU,-1)+speye(size(LU));
U=triu(LU); 

if (s == 1)
    sA = R\A;
else
    sA = A;
end

myErr = lu_normest(sA(rowp,colp),L,U)/norm(A,1)
%umfErr = lu_normest(D(:,p)\A(:,q),l,u)
umfErr = lu_normest(p*(r\A)*q,l,u)/norm(A,1)

umfpnnz = nnz(l)+nnz(u) - size(A,1)
mynnz = nnz(LU) %+ nnz(paddingZ)

myflop = luflop(L,U)


umfflop = luflop(l,u)

flopratio = myflop/umfflop

hardware_flp_cnt 

ratio = hardware_flp_cnt / myflop

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
