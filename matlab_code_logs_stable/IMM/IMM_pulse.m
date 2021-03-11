function track = IMM_pulse(z,track,par)
%      par - struct containing track parameters 
%      track - struct containing tracking variables 
%    Output parameter is:
%      track - the updated track
warning off
A_m=par.A_m; %transition prob matrix
F=par.F; % mode dinamics

%combined estimates for output
track.x=zeros(2,1);
track.P=zeros(2);
M=par.models;
% state propagation
mu_hat=track.mu*A_m;      
track.mu_hat=mu_hat;
B_m=A_m.*repmat(track.mu',1,M)./repmat(mu_hat,M,1);%  backwards transition probability matrix
B_m(isnan(B_m)|isinf(B_m))=0;
Lambda=zeros(1,M);
zi=zeros(1,M);
S=zeros(1,M);
Pi=zeros(2,2,M);
for k=1:M
	xm=track.xm*B_m(:,k);
	%mix  xm and Pm
	Pm=-xm*xm';
	for n=1:M
		Pm_m=track.Pm(:,:,n)+track.xm(:,n)*track.xm(:,n)';
		Pm=Pm+B_m(n,k)*Pm_m;
	end
	track.xm(:,k) = F(:,:,k)*xm;
	Pi(:,:,k)= F(:,:,k) * Pm * F(:,:,k)' + par.Q(:,:,k); %
	S(k) = par.H*Pi(:,:,k)*par.H' + par.R(k); %
	% zhat is the predicted measurement
	zhat=par.H*track.xm(:,k);    %
	zi(k)=z-zhat;%innovation
	Lambda(k)=max(eps,1/sqrt(2*pi)/sqrt(S(k))*exp(-1/2*zi(k)^2/S(k)));
end

track.zi=zi;
track.S=S;
Pz=Lambda*mu_hat';
for k=1:M
	track.mu(k)=Lambda(k)*mu_hat(k)/Pz;
	%kalman gain
	K=Pi(:,:,k)*par.H'/S(k);
	%state estimate update
	track.xm(:,k)=track.xm(:,k)+K*zi(:,k);
	%error covariance matrix update
	track.Pm(:,:,k)=Pi(:,:,k)-K*S(k)*K';
	track.x=track.x+track.mu(k)*track.xm(:,k);
	track.P=track.P+track.mu(k)*(track.Pm(:,:,k)+track.xm(:,k)*track.xm(:,k)');
end
track.P=track.P-track.x*track.x';
warning on