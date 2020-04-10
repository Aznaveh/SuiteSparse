%T = load ('../Demo/sinc18');
%T = load ('../Demo/sinc15.out');
%T = load ('../Demo/umfSinc15.out');
%T = load ('../Demo/umfSinc18');
%T = load ('../Demo/twoton');
T = load ('../Demo/umftwoton');
M = T(:,1);
N = T(:,2);
K = T(:,3);

[KK, ia, ic] = unique(K,'sorted');
F = zeros(size(KK));
for  i = 1:length(K)
    ind = ic (i);
    F(ind) = F(ind)+ M(i)*N(i)*K(i);
end

[NN, nia, nic] = unique(N,'sorted');
NF = zeros(size(NN));
for  i = 1:length(N)
    ind = nic (i);
    NF(ind) = NF(ind)+ M(i)*N(i)*K(i);
end



bar(KK, F, 4)
%bar(NN, NF, 10)
