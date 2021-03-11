function [F,mu,v,f0,mu0,f,L]=iterVM(z,fran)
%iterative algorithm for computing the frequency from irregular times of
%arrival
% d=4.66124e-13;
mu=0;
N=2^16;
f=linspace(fran(1),fran(2),N+1);
F=mean(fran);
clear mu0 f0
for k=1:20
	mu0(k)=mu;
	f0(k)=F;
	L=sum(cos(2*pi*(z*f-mu))); %VM periodogram
	[v,k0]=max(L);
	F=f(k0);
	y=exp(2*pi*1i*z*F);
	mu=angle(mean(y))/2/pi;
% 	if mu==mu0(k) & F==f0(k)
% 		break
% 	end
	% 	kappa=1/(1-abs(mean(y)));
	f=F+linspace(-1/2^(k+1),1/2^(k+1),N+1);
end
