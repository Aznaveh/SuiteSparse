function drawFronts(A)
e=1;
w=2;

[m n]= size(A)

mn= min(m,n)

%pcol=colamd(A);
pcol=metis(A,'col');
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


p=etree(A,'col');

etreeplot(A'*A);
hold on

pause;
disp('Blue is zero part, Red is nonzero box. Press a key to Continue...');

%find first of j
f=zeros(n,1)-1;
for i=1:n
    if p(i)== 0 %the etree is either forest or we are in the root
        if f(i) == -1 
            f(i)= i;
        end
    elseif f(i)== -1
        f(i)=i;
        t=p(i);
        while t > 0 && f(t) == -1     %traversing up the tree
            f(t) = i;
            t=p(t);
        end
    end
end

%findign fstElinCl and lstElinCl
lstElinCl=ones(n,1);
fstElinCl=(m+1)*ones(n,1);
[I,J,~]=find(A);
for k=1: nnz(A)
    i= I(k); j=J(k);
    if j <= m
        fstElinCl(j) = min(fstElinCl(j),i);
        lstElinCl(j) = max(lstElinCl(j),i);
    end
end


for j = 1:mn
    spy(A);
    hold on
    drawbox(1,fstElinCl(f(j)),f(j),j,'b',w,e);
    drawbox(fstElinCl(f(j)),lstElinCl(j),1,f(j),'b',w,e);
    drawbox(fstElinCl(f(j)),lstElinCl(j),f(j),j,'r',w,e);
    pause(min(.5,500/(m*n) ));
    hold off
end
