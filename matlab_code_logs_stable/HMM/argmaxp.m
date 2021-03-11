function j=argmaxp(L)
%
N=size(L,1);
n=(N+1)/2;
L(1:n-1,:)=0;
[~,j]=nanmax(L);