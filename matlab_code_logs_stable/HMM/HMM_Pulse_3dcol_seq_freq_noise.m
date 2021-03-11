function [alpha,evidence] = HMM_Pulse_3dcol_seq(z,fran,fdran,glitch,alpha_m0,evidence_m0,s,M)
% Forward-backward --- used when are searching for glitches
% keep alpha(0:n-1) and beta(n:T) computed with M_0 model
% keep evidence(0:n-1) 
%   
% Inputs:
% z- TOAs
% fran - frequency bins
% fdran -fdot bins  
% glitch - glitch position
% s -likelihood
% Outputs:
% alpha - forward probability
% evidence - evidence at each time step

%T -time chunks, Q- fdot bins; V-f bins
% tic
global sigma 
T=length(z);
K=length(fdran);
N=length(fran);
% sigma = 1.0e-16; %process noise
df=fran(2)-fran(1);
dfdot=fdran(2)-fdran(1);
k=(K-1)/2;
fdran_ext=dfdot*(-2*k:2*k)+fdran((K+1)/2);% extend f_dot
%allocate large variables
alpha=alpha_m0;
evidence=evidence_m0(1:M);
%forward probabilities
alpha0=alpha_m0(:,:,M);
alpha(:,:,1:M)=alpha_m0(:,:,1:M);
for n=M:T
	gi=find(any(n==glitch));
	[q,r,mu]=fokker_plank_f_freq_noise(sigma,z(n),fran,fdran_ext,df,dfdot,gi);
	[y]=colmaxf2d_r(alpha0,q,mu,r,gi);
	alpha0=y + s(:,:,n);%bayes rule
	evidence(n+1)= logsumexp(alpha0(:));
	alpha(:,:,n+1)=alpha0;
end
% toc
