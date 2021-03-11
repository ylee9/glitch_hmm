function [path,delta,mode,psi,phi,alpha,evidence,pg, residuals] = viterbi_Pulse_3dcol(z,f0,fdot,glitch)
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
K=length(fdot);
N=length(f0);

% delta = zeros(Q,V);
psi = zeros(T,K,N,2);
phi=zeros(T,K,N);
path = zeros(T,2); % 2D path (fdot and f0)
mode=zeros(1,T);
pmode=zeros(K,N,2,T);
residuals = zeros(1, T);
global kappa
lik=cos(2*pi*(z(1)*f0+1/2*z(1).^2*fdot'));
log_besseli = log(besseli(0, kappa));
if(log_besseli == Inf)
    % For large kappa, I_0(kappa) ~ exp(kappa)/sqrt(2pi*kappa)
    log_besseli = kappa - 0.5*log(2*pi*kappa);
end

% lik=WGlik(f0,fdot,z(1),v);
delta = lik;
alpha=kappa*lik - log_besseli - log(2*pi);
evidence(1)=logsumexp(alpha(:));
alpha=alpha - evidence(1);
% disp('forward')
global sigma;
df=mean(diff(f0));
dfdot=mean(diff(fdot));
fr=fdot(end)-fdot(1);% range of fdot
%fdot_ext=fdot(1)-fr/2-dfdot:dfdot:fdot(end)+fr/2+dfdot;% extend f_dot
k = (K-1)/2;
fdot_ext=dfdot*(-2*k:2*k)+fdot((K+1)/2);
for t=2:T
    gi=find(any(t==glitch));
    [q, r, mu] = fokker_plank_f(sigma, z(t), f0, fdot_ext, df, dfdot, gi);

	lik=cos(2*pi*(z(t)*f0+1/2*z(t).^2*fdot'));

    [y,mc,modd,alpha0,pmode]=colmaxf2d(delta,q,mu,alpha,gi,r);
	%bayes rule
	delta=y+lik;%	
	pmode(:,:,1)=pmode(:,:,1)+lik;
	pmode(:,:,2)=pmode(:,:,2)+lik;
    pg_all(:,:,t) = pmode(:, :, 2) - logsumexp(pmode, 3);
	
	log_s = kappa*lik - log_besseli - log(2*pi);
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
	try
	path(t,:) = squeeze(psi(t+1,path(t+1,1),path(t+1,2),:,:));
	mode(t)   = squeeze(phi(t+1,path(t+1,1),path(t+1,2)));
    phase = (z(t)*f0(path(t, 2))+1/2*z(t).^2*fdot(path(t, 1)));
    residuals(t) = phase - round(phase);
	pg(t)=squeeze(pg_all(path(t+1,1),path(t+1,2),t+1));
	catch
	end
end
%  toc