function k=modm(k,n)
% modulo without zero
k=mod(k,n);
k(k==0)=n;