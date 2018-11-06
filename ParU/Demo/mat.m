%%%%%%%% START %%%%%%%%%%%%%%%%%%% 
clear all
% A  is  4 x 4 
LU = zeros(4,4);
npivots =[]; 
S = zeros(4,4);
% nf=2
% anz = 9  rjsize=5
S(1,1)= 4.0000000000000000;
S(1,4)= -1.0000000000000000;
S(2,1)= 1.0000000000000000;
S(2,4)= 4.0000000000000000;
S(3,2)= 4.0000000000000000;
S(3,3)= 1.0000000000000000;
S(3,4)= 1.0000000000000000;
S(4,2)= -1.0000000000000000;
S(4,3)= 4.0000000000000000;
%~~~~~~~  Assemble Front 0 start ~~~~
% fp=1 pivotal columns:clo1=0...col2=0
%L part:
cols{1} = [1 ];
rows{1} = [1 2 ];
Luf{1}= [  4.0000000000000000 ;
     0.2500000000000000 ;
   ];
Us{1} =[];
Ucols{1}=[];
Urows{1}=[];
rowCount=2;
% U part After TRSM: 1 x 1
Ucols{1} = [4 ];
Urows{1} = [1 ];
Us{1} = [ -1.0000000000000000 ;
    ];
%~~~~~~~Assemble Front 0 finished
%~~~~~~~  Assemble Front 1 start ~~~~
% fp=3 pivotal columns:clo1=1...col2=3
%L part:
cols{2} = [2 3 4 ];
rows{2} = [3 4 2 ];
Luf{2}= [  4.0000000000000000  1.0000000000000000  1.0000000000000000 ;
     -0.2500000000000000  4.2500000000000000  4.2500000000000000 ;
     0.0000000000000000  0.0000000000000000  0.2500000000000000 ;
   ];
Us{2} =[];
Ucols{2}=[];
Urows{2}=[];
err = 1e-12; [m n] =size(S);
oldR=[]; c=[];
%Finalizing the permutation
for f=1:2
	npivots(f) = length(cols{f});
	oldR=[oldR, rows{f}(1:npivots(f))]; c=[c, cols{f}];
end
oldR = [oldR setdiff(1:m,oldR)];
newR(oldR)=1:length(oldR);
for f=1:2
	LU(newR(rows{f}),cols{f})=Luf{f};
	LU(newR(Urows{f}),Ucols{f})=Us{f};
end
L=tril(LU,-1)+eye(size(LU));
U=triu(LU); U=U(1:n,:);
spparms('spumoni',3);
fprintf('Matlab\n');
[l,u,p]=lu(S);
norm(p*S-l*u)
fprintf('Paru\n');
norm(S(oldR,c)-L*U)
if( (norm(S(oldR,c)-L*U)) < err )
	fprintf('Pass\n')
else
	fprintf('Fail\n')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*
