
matrix_size = 100;
band_width = 3;
k = 1; % 1/k part of the first row is nnz

A = spdiags(ones(matrix_size,2*band_width+1),...
            -band_width:band_width,matrix_size,matrix_size);
A = A + 4*speye(matrix_size);
%A (1,:) = 10*ones(1,matrix_size);
%A (end,:) = 10*ones(1,matrix_size);
A (1,1:matrix_size/k) = 10*ones(1,matrix_size/k);
mmwrite('../Matrix/ParUTst/tmp.mtx', A);
