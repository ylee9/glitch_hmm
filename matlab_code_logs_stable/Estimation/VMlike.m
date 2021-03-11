function [L,J,H]=VMlike(par,x)
%
% z=[0;diff(x)];
w=ones(size(x))'/length(x);
f=par(1);df=par(2);mu=par(3);
% ddf=par(4);
ddf=0;
L=-w*cos(2*pi*(x.*f-x.^2*df/2-mu));
J=w*[ -2*x*pi.*sin(2*pi*(mu - x.*(f - (df*x)/2))), x.^2*pi.*sin(2*pi*(mu - x.*(f - (df*x)/2))), 2*pi*sin(2*pi*(mu - x.*(f - (df*x)/2)))];
H=[  w*(4*x.^2*pi^2.*cos(2*pi*(mu - x.*(f - (df*x)/2)))), -w*(2*x.^3*pi^2.*cos(2*pi*(mu - x.*(f - (df*x)/2)))),  -w*(4*x*pi^2.*cos(2*pi*(mu - x.*(f - (df*x)/2))));...    
 -w*(2*x.^3*pi^2.*cos(2*pi*(mu - x.*(f - (df*x)/2)))),    w*(x.^4*pi^2.*cos(2*pi*(mu - x.*(f - (df*x)/2)))), w*(2*x.^2*pi^2.*cos(2*pi*(mu - x.*(f - (df*x)/2))));...
 -w*(4*x*pi^2.*cos(2*pi*(mu - x.*(f - (df*x)/2)))),  w*(2*x.^2*pi^2.*cos(2*pi*(mu - x.*(f - (df*x)/2)))),     w*(4*pi^2*cos(2*pi*(mu - x.*(f - (df*x)/2))))];
