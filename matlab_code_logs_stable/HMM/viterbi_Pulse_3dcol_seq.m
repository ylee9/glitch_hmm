function [path,delta,mode,psi,phi,evidence,delta0,evidence0] = viterbi_Pulse_3dcol_seq(z,f0,fdot,glitch,delta,evidence)
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
 tic
%T -time chunks, Q- f dot bins; V-f0 bins
T=length(z);
Q=length(fdot);
V=length(f0);

% delta = zeros(Q,V);
psi = zeros(T,Q,V,2);
phi=zeros(T,Q,V);
path = zeros(T,2); % 2D path (fdot and f0)
mode=zeros(1,T);
% v=2*pi*f*2e-3;%
% lik=WGlik(f0,fdot,z(1),v);
% delta = lik;
kappa=10;
scale=1;
% evidence=exp(kappa*lik)/2/pi/besseli(0,kappa)/scale;
% disp('forward')
sigma = 1.0e-16; %process noise
df=mean(diff(f0));
dfdot=mean(diff(fdot));
% pmode=zeros(Q,V,2,T);
% pmode(:,:,1,1)=1;
for t=1:T
% 	d=squeeze(delta(t-1,:,:));
	m=max(1,min(V,ceil(sigma*z(t).^(3/2)/3^(1/2)*(2*6+1)/df)));
	if (m/2)==round(m/2), m=m-1;end
	[q]=fokker_plank(sigma,z(t),f0(1:m),fdot,df,dfdot);
	mu=[f0-z(t)*fdot']; %map to the next level
	mun=round((mu-f0(1))/df)+1; % 
 	gi=any(glitch==t);
	if gi
		delta0=delta;
		evidence0=evidence;
	end
	[y,mc,modd,evidence,pmode(:,:,:,t)]=colmaxf2d(delta,q,mun,evidence,gi);
	%bayes rule
	lik=cos(2*pi*(z(t)*f0+1/2*z(t).^2*fdot'));

%     lik=WGlik(f0,fdot,z(t),v);
% 	delta(t,:,:)=y+lik;
	delta=y+lik;%-min(lik(:));
	evidence=(evidence).*exp(kappa*lik)/2/pi/besseli(0,kappa)/scale;
% 	disp(sum(evidence(:)))
% 	disp(t)
	psi(t,:,:,:)=mc;
	phi(t,:,:)=modd;
%  	disp(t)
end

% disp('backward')
% tmp=squeeze(delta(T,:,:));
S=max(delta(:));

[path(T,1),path(T,2)]=find(delta==S,1,'first');
[mode(T)] = phi(T,path(T,1),path(T,2)); 
for t=T-1:-1:1
	try
	path(t,:) = squeeze(psi(t+1,path(t+1,1),path(t+1,2),:,:));
	mode(t)   = squeeze(phi(t+1,path(t+1,1),path(t+1,2)));
	catch
	end
end
 toc