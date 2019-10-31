function save_sp_mat(A)
%    [row col v] = find(A);
%    T = [size(A) nnz(A)];
%    T = [T; row col v];
%    dlmwrite('../Matrix/ParUTst/tmp.mtx',T,'delimiter', '\t')
mmwrite('../Matrix/ParUTst/tmp.mtx',A)

