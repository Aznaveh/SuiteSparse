[m n]= size(A)
mn= min(m,n);

pcol=metis(A,'col');
%pcol=colamd(A);
A=A(:,pcol);
[p,q]=etree(A,'col');


%making staircase
leftmost=(n+1)*ones(n,1);
[I,J,~]=find(A);
for k=1: nnz(A)
    i= I(k); j=J(k);
    if i <= n
        leftmost(i) = min(leftmost(i),j);
    end
end
[sortedleftm,prow]=sort(leftmost);
A=A(prow(1:mn),:);
