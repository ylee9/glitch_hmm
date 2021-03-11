function [L,P]=WGlik(f0,fdot,z,v)	
phi=mod(z*f0+1/2*z.^2*fdot',1);
P=zeros(size(phi));
for n=-100:100
	P=P+1/sqrt(2*pi)/v*exp(-2*pi^2/v^2*(phi+n).^2);
end
L=log(P);
	
