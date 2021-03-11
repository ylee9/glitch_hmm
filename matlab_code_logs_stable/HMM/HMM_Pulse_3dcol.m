function [path,gamma,alpha,beta,evidence,s, residuals] = HMM_Pulse_3dcol(z,fran,fdran,glitch)
% Forward-backward 
%   
% Inputs: o
% z- TOAs
% fran - frequency bins
% fdran -fdot bins  
% glitch - glitch position
% Outputs:
% alpha - forward probability
% evidence - evidence at each time step

%T -time chunks, Q- fdot bins; V-f bins
global sigma
global kappa_per_toa
%sigma
T=length(z);
K=length(fdran);
length(fdran);
N=length(fran);
% sigma = 1.0e-16; %process noise
df=fran(2)-fran(1);
dfdot=fdran(2)-fdran(1);
k=(K-1)/2;
fdran_ext=dfdot*(-2*k:2*k)+fdran((K+1)/2);% extend f_dot
%allocate large variables
alpha=zeros(K,N,T+1);
beta=zeros(K,N,T+1);
s=zeros(K,N,T);
residuals = zeros(1, T);
% gamma=zeros(K,N,T);
path=zeros(T,2);
evidence=zeros(1,T+1);
%forward probabilities
alpha0 = zeros(K,N);
evidence(1)=logsumexp(alpha0(:));
%alpha0=alpha0 - evidence(1);
alpha(:,:,1)=alpha0;
% draw=0;

global f_fiducial
global fd_fiducial
global fdd_fiducial;
f_fiducial_running = f_fiducial;
fd_fiducial_running = fd_fiducial;
t = zeros(T,1);

for n=1:T
    kappa = kappa_per_toa(n);
    gi=find(any(n==glitch));
    [q,r,mu]=fokker_plank_f(sigma,z(n),fran,fdran_ext,df,dfdot,gi);
    [y]=colmaxf2d_r(alpha0,q,mu,r,gi);
    f_fiducial_running = f_fiducial_running - fd_fiducial_running*z(n) + 0.5*fdd_fiducial*z(n)^2;;
    fd_fiducial_running = fd_fiducial_running - fdd_fiducial*z(n);
    x_fiducial = 2*pi*(z(n)*f_fiducial_running + 0.5*fd_fiducial_running*z(n)^2 + 1/6*fdd_fiducial*z(n)^3);

     x = 2*pi*(z(n)*fran+1/2*z(n)^2*fdran') + x_fiducial;
    %lik=cos(2*pi*(z(n)*fran+1/2*z(n)^2*fdran'));
    log_besseli = log(besseli(0, kappa));
    if(log_besseli == Inf)
        % For large kappa, I_0(kappa) ~ exp(kappa)/sqrt(2pi*kappa)
        log_besseli = kappa - 0.5*log(2*pi*kappa);
    end
    log_s = kappa*cos(x) - log_besseli - log(2*pi);
%    log_s = log_s - logsumexp(log_s(:));
    s(:,:,n)= log_s;
    alpha0=y + s(:,:,n);%Bayes rule
    evidence(n+1) = logsumexp(alpha0(:));
    alpha(:,:,n+1)=alpha0;
end

% backward probabilities
beta0=zeros(K,N);
beta0=beta0 - evidence(T+1);
beta(:,:,T+1)=beta0;%T-th beta
for n=T:-1:2
	gi=find(any(n==glitch));
	y=beta0 + s(:,:,n);
	[q,r,mu]=fokker_plank_b(sigma,z(n),fran,fdran_ext,df,dfdot,gi);
	[beta0]=colmaxf2d_r(y,q,mu,r,gi);

	beta(:,:,n)=beta0;
end
gamma=alpha(:,:,2:T+1) + beta(:,:,2:T+1);
f_fiducial_running = f_fiducial;
fd_fiducial_running = fd_fiducial;
for n=1:T
%   if draw
%       pcolor(fran,fdran,gamma(:,:,n)),shading interp
%       title(num2str(n))
%       drawnow
%   end
    gamma(:,:,n) = gamma(:,:,n) - logsumexp(logsumexp(gamma(:,:,n)));
    [path(n,1),path(n,2)]=find(gamma(:,:,n)==max(max(gamma(:,:,n))),1);
    f_fiducial_running = f_fiducial_running - fd_fiducial_running*z(n) + 0.5*fdd_fiducial*z(n)^2;
    fd_fiducial_running = fd_fiducial_running - fdd_fiducial*z(n);
    phase_fiducial = (z(n)*f_fiducial_running + 0.5*fd_fiducial_running*z(n)^2 + 1/6*fdd_fiducial*z(n)^3);


    phase = (z(n)*fran(path(n, 2))+1/2*z(n).^2*fdran(path(n, 1))) + phase_fiducial;
    residuals(n) = phase - round(phase);
end
