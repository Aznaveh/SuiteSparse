function x = lsolve(L,b)
%return L\b
[m n] = size (L); % assert(m==n)
x = b;
step = 3;
k=1:step:n;

k=1:step:n;
if (k(end) ~= n)
    k = [k,n]
end

for j = 2:size(k,2)
    j1 = k(j-1)+1;
    if (j1==2) 
        j1=1;
    end
    j2 = k(j);


    x(j1:j2) = L(j1:j2,j1:j2)\x(j1:j2); %trsv
    x(j2+1:n) = x(j2+1:n) - L(j2+1:n,j1:j2)*x(j1:j2); %dgemv
end


