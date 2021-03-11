function [y]=colmaxf2d_r(alpha,Q,mu,R,gi)
% Q = fokker_plank(sigma,z(t),f0(1:m(t)),fdot); t - measurement
% number,m=min(3001,ceil(sigma*z(t)^(3/2)/3^(1/2)*(2*n+1)/df)); for n=6 standard
% diviations, upper limit for 11 fdot bins is ~ 3001 ----
% mun=max(0,round((mu-f0(1))/df)+1) and mu=[f0-z(t)*fdot'];  
% L -glitch size (in freq bins)

[K,N]=size(alpha);
%dynamics for Mode A
M=size(Q,2);
y0=[ones(K,(M-1)/2)*-Inf alpha ones(K,(M-1)/2)*-Inf];
al=im2col(y0,[K M],'sliding');
alphaL=ones(K,N)*-Inf;
alphaG=ones(K,N)*-Inf;
% e=0.0;
% A=[1-e e];
% if gi %glitch indicator
% 	A=[e 1-e];
% end 
%dynamics for Mode L
if isempty(gi)
for k=1:K % loop over predicted fdot
	Qn=shift_q(Q,k);
	P=repmat(Qn(:),1, N);
	D=logsumexp(al + P);
	%work out the correct displacement and write it into arrays
	place=mu(k)+(1:N);
	i=find(place>0 & place<N+1);
	alphaL(k,place(i))=D(i);
end
y=alphaL;
else
%dynamics for Mode G
for n=1:N
    
	Rn=R(:,(1:N)+N-n);
  	%Rn=Rn/sum(Rn(:));
	loglikes = alpha + Rn;
    alphaG(:,n) = logsumexp(loglikes(:));
end
y=alphaG;
end
%combine the scores
% y=A(1)*alphaL+A(2)*alphaG;

