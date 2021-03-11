function [theta,n,theta0]=em_pulse_plus(dz,theta)
% P=0.60688662425;
% Pdot=4.66124e-13*P;
N=length(dz);
% theta=[P;Pdot];
% assume geometric distribution for n with parameter 1/lambda, also works for poisson 
for k=1:10000
	theta0(:,k)=theta;
	p=theta(1);
	dp=theta(2);
	%solve quadratic equation for n
	n=round(-(dp+ 2*p - 2*(dp^2/4 + dp*p + 2*dz*dp + p^2).^(1/2))/(2*dp));
	A=[sum(n.^2) sum((n+1).*n.^2)/2; sum((n+1).*n.^2)/2 sum((n+1).^2.*n.^2)/4];
	b=[sum(dz.*n); sum(dz.*(n+1).*n/2)];
	theta=A\b;
	if norm(theta0(:,k)-theta)==0
		disp(k)
		break
	end
end