function [theta,n,theta0]=em_pulse(dz,P)
% P=0.60688662425;
% Pdot=4.66124e-13*P;
N=length(dz);
Z=sum(dz)/N;
theta=[P];
% assume geometric distribution for n with parameter 1/lambda, also works for poisson 
for k=1:10000
	theta0(k)=theta;
	p=theta(1);
    n=round(dz/p);
	m1=sum(dz.*n)/N;
	m2=sum(n.^2)/N;
	lambda=sum(n)/N;
	p=Z/lambda;
	theta=[p];
	if theta0(k)==theta
		disp(k)
		break
	end
end