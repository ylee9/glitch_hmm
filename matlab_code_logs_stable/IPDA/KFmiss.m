function [Y,track]=KFmiss(meas,z)
% kalman filter with missing data and unknown label
% p=0.60688662425;
% dp=4.6612e-13*p;
P=eye(2);
% R=5e-6;
H=[1 0];
q=1e-32;

track=struct('x',[],'xhat',[],'P',[],'Phat',[]);
track_tmp=struct('x',[],'xhat',[],'P',[],'Phat',[]);

for n=1:length(meas{1})
	track(n).x=[meas{1}(n);-1.25e-12];
	track(n).P=P;
end
Y(:,:,1)=[track.x];
for ii=2:length(z)
	[Fn,Qn]=comp_noise_cov(z(ii),q);
	R=1e-6;
	k=1;
	for m=1:length(track)
		track(m).xhat=Fn*track(m).x;
		track(m).Phat=Fn*track(m).P*Fn'+Qn;
		S  = H*track(m).Phat*H' + R;
		%data association
		zi = meas{ii}-H*track(m).xhat;
		[~,n]=min(zi.^2/S);
		zi=zi(n);
		K = track(m).Phat*H'/S;
        P=track(m).Phat-K*H*track(m).Phat;
		for n=1:length(zi)
			%kalman gain
			track_tmp(k).x = track(m).xhat+K*zi(n);
			track_tmp(k).xhat=track(m).xhat;
			%error covariance matrix update
			track_tmp(k).P = P;
			track_tmp(k).Phat=track(m).Phat;
			k=k+1;
		end
	end
	track=track_tmp;
	for n=1:length(track)
	Y(:,n,ii)=track(n).x;
	end
end
	

