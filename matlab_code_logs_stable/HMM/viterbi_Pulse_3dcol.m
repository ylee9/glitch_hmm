function [path,delta,mode,psi,phi,alpha,evidence,pg, residuals] = viterbi_Pulse_3dcol(z,fran,fdran,glitch)
% VITERBI Find the most-probable (Viterbi) path through the HMM state trellis.
% path = viterbi(z,f0,fdot)
%   
% Inputs:
% 
%
% Outputs:
% path(t) = q(t), where q1 ... qT is the argmax of the above expression.
% delta(j,t) = prob. of the best sequence of length t-1 and then going to state j, and O(1:t)
% psi(j,t) = the best predecessor state, given that we ended up in state j at t
%  tic
%T -time chunks, Q- f dot bins; V-f0 bins
T=length(z);
K=length(fdran);
N=length(fran);

% delta = zeros(Q,V);
psi = zeros(T,K,N,2);
phi=zeros(T,K,N);
path = zeros(T,2); % 2D path (fdot and f0)
mode=zeros(1,T);
pmode=zeros(K,N,2,T);
residuals = zeros(1, T);
global kappa_per_toa
kappa = kappa_per_toa(1)
log_besseli = log(besseli(0, kappa));
if(log_besseli == Inf)
    % For large kappa, I_0(kappa) ~ exp(kappa)/sqrt(2pi*kappa)
    log_besseli = kappa - 0.5*log(2*pi*kappa);
end
lik=kappa*cos(2*pi*(z(1)*fran+1/2*z(1).^2*fdran'));


% lik=WGlik(f0,fdot,z(1),v);
delta = lik;
alpha=lik;
evidence(1)=logsumexp(alpha(:));
alpha=alpha - evidence(1);
% disp('forward')
global sigma;
df=mean(diff(fran));
dfdot=mean(diff(fdran));
fr=fdran(end)-fdran(1);% range of fdot
%fdot_ext=fdot(1)-fr/2-dfdot:dfdot:fdot(end)+fr/2+dfdot;% extend f_dot
k = (K-1)/2;
fdot_ext=dfdot*(-2*k:2*k)+fdran((K+1)/2);

global f_fiducial
global fd_fiducial
global fdd_fiducial;
fd_fiducial_running = fd_fiducial - fdd_fiducial*z(1);
f_fiducial_running = f_fiducial - fd_fiducial_running*z(1)

for t=2:T
    kappa = kappa_per_toa(t)
    gi=find(any(t==glitch));
    [q, r, mu] = fokker_plank_f(sigma, z(t), fran, fdot_ext, df, dfdot, gi);

    f_fiducial_running = f_fiducial_running - fd_fiducial_running*z(t) + 0.5*fdd_fiducial*z(t)^2;;
    fd_fiducial_running = fd_fiducial_running - fdd_fiducial*z(t);
    x_fiducial = 2*pi*(z(t)*f_fiducial_running + 0.5*fd_fiducial_running*z(t)^2 + 1/6*fdd_fiducial*z(t)^3);
    x = 2*pi*(z(t)*fran+1/2*z(t)^2*fdran') + x_fiducial;
    
	lik=kappa*cos(x);
    [y,mc,modd,alpha0,pmode]=colmaxf2d(delta,q,mu,alpha,gi,r);
	%bayes rule
	delta=y+lik;%	
	pmode(:,:,1)=pmode(:,:,1)+lik;
	pmode(:,:,2)=pmode(:,:,2)+lik;
    pg_all(:,:,t) = pmode(:, :, 2) - logsumexp(pmode, 3);
	
	log_s = lik;
	alpha=alpha0 + log_s;
    evidence(t) = logsumexp(alpha(:));
	alpha=alpha - evidence(t);
	psi(t,:,:,:)=mc;
	phi(t,:,:)=modd;
end

S=max(delta(:));
[path(T,1),path(T,2)]=find(delta==S,1,'first');
[mode(T)] = phi(T,path(T,1),path(T,2)); 
for t=T-1:-1:1
	
	path(t,:) = squeeze(psi(t+1,path(t+1,1),path(t+1,2),:,:));
	mode(t)   = squeeze(phi(t+1,path(t+1,1),path(t+1,2)));
    
	pg(t)=squeeze(pg_all(path(t+1,1),path(t+1,2),t+1));
end

fd_fiducial_running = fd_fiducial;
f_fiducial_running = f_fiducial;
for n = 1:T
    f_fiducial_running = f_fiducial_running - fd_fiducial_running*z(n) + 0.5*fdd_fiducial*z(n)^2;
    fd_fiducial_running = fd_fiducial_running - fdd_fiducial*z(n);
    phase_fiducial = (z(n)*f_fiducial_running + 0.5*fd_fiducial_running*z(n)^2 + 1/6*fdd_fiducial*z(n)^3);

    phase = (z(n)*fran(path(n, 2))+1/2*z(n).^2*fdran(path(n, 1))) + phase_fiducial;
    residuals(n) = phase - round(phase);
end
%  toc
