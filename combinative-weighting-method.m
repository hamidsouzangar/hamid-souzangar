clc
clear all
%%
% The codes of  Combinative weighting method
% Hamidreza Souzangarzadeh
%Multi-objective optimization of cylindrical segmented tubes as energy 
%absorbers under oblique crushes: D-optimal design and integration of
%MULTIMOORA with combinative weighting
try
	% See if there is an existing instance of Excel running.
	% If Excel is NOT running, this will throw an error and send us to the catch block below.
	Excel = actxGetRunningServer('Excel.Application');
	% If there was no error, then we were able to connect to it.
	Excel.Quit; % Shut down Excel.
catch
	% No instance of Excel is currently running.
end
input=xlsread('weight','input');% NUmber of each tube

m=16
n=10
for j=1:n
ii=zeros(m,n);
for i=1:m;
  i1(i,j)= input(i,j)-input(17,j); % CHECK CHECK CHECK
  % i1(i,2)=  input(i,2)-input(9,2);

end 
end
i2=i1.*i1;
ii=zeros(n,n);
ii2=zeros(n,n);
mak=zeros(n,n)
for j=1:n;
    for jj=1:n;
    for i = 1:m;
      if j==jj;
o=1 
else 
o=-1
end
      ii (j,jj,i)= i1(i,j)*i1(i,jj);
      
    end
    ii2_1(jj)=sum(i2(:,jj))';
    sor(j,jj)=sum(ii(j,jj,:))*input(18,j)*input(18,jj); % CHECK CHECK CHECK
        mak(j,jj)=sqrt(ii2_1(j)*ii2_1(jj));

    end
end
   
   R=sor./mak
   %%
   
   yek = ones(n,n);
   R2 = yek-R;
   sor2 = sum(R2(:,1:n));
   mak2 = sum(sor2);
   Wc = sor2./mak2;
   
   %%
   input2 = xlsread('weight','input','a2:j17')
   p= sum (input2(1:m, :));
   P2 = repmat(p,m,1);
      Pij = input2./P2;
      E=-(sum(Pij(1:m,:).*log(Pij)))/log(m);
      E1=ones(1,n)-E;
      Wo=(E1)./sum(E1);
      %%
      Ws = [0.1033333333333330 0.1033333333333330 0.1368888888888890 0.0973333333333333 0.0877777777777778 0.1282222222222220 0.0877777777777778 0.0711111111111111 0.1202222222222220 0.0640000000000000 
         ];
W1 = (Wc.*Wo.*Ws).^(1/3);
Wt=W1/sum(W1)
