function x = usolve(U,b)
%return U\b
[m n] = size (U); % assert(m==n)
x = b;
step = 3;

k=1:step:n;
if (k(end) ~= n)
    k = [k,n];
end

for j = size(k,2):-1:2
    j1 = k(j-1)+1;
    if (j1==2) 
        j1=1;
    end
    j2 = k(j);


    x(j1:j2) = U(j1:j2,j1:j2)\x(j1:j2); %trsv
    x(1:j1-1) = x(1:j1-1) - U(1:j1-1,j1:j2)*x(j1:j2); %dgemv
end


