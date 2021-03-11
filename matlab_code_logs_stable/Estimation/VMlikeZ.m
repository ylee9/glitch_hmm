function [L,J,H]=VMlikeZ(par,x)
%
z=[0;diff(x)];
w=ones(size(z))'/length(z);
f=par(1);df=par(2);mu=par(3); ddf=par(4);
L=-w*cos(2*pi*(z.*f-z.*x.*df-mu-z.*x.^2/2*ddf));
J=[ -2*w*(z*pi.*sin(2*pi*(mu - f*z + df*z.*x))), 2*w*(z.*x*pi.*sin(2*pi*(mu - f*z + df*z.*x))), 2*w*(pi*sin(2*pi*(mu - f*z + df*z.*x)))];
H=[    4*w*(z.^2*pi^2.*cos(2*pi*(mu - f*z + df*z.*x))),  -4*w*(z.^2.*x*pi^2.*cos(2*pi*(mu - f*z + df*z.*x))),  -4*w*(z*pi^2.*cos(2*pi*(mu - f*z + df*z.*x)));
 -4*w*(z.^2.*x*pi^2.*cos(2*pi*(mu - f*z + df*z.*x))), 4*w*(z.^2.*x.^2*pi^2.*cos(2*pi*(mu - f*z + df*z.*x))), 4*w*(z.*x*pi^2.*cos(2*pi*(mu - f*z + df*z.*x)));
     -4*w*(z*pi^2.*cos(2*pi*(mu - f*z + df*z.*x))),     4*w*(z.*x*pi^2.*cos(2*pi*(mu - f*z + df*z.*x))),     4*w*(pi^2*cos(2*pi*(mu - f*z + df*z.*x)))];

