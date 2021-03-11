function tracker = LIPDAF(z_old,z,tracker_par,tracker,tracker_init_val)
% LIPDAF Linear Integrated Probability Data Association Filter
%    TRACKER = IPDAF(WFMS,Z_OLD,Z,TRACKER_PAR,TRACKER,CLUTTERMAPPER_PAR,TRACKER_INIT_VAL)
%    This is the IPDAF as described in the report.
%    Input parameters are:
%      Z_OLD - the old measurements
%      Z - the measurements as seen by the radar
%      TRACKER_PAR - struct containing tracker parameters (see INIT_PARS)
%      TRACKER - struct containing tracking variables (see INIT_VARS)
%      TRACKER_INIT_VAL - struct containing the values of the tracker
%                         for initialization purpose
%    Output parameter is:
%      TRACKER - the updated TRACKER

warning off
%----------------------
% initialize and augment tracks
tracker=init_tracks(z,z_old,tracker_par,tracker,tracker_init_val);
% kill similar tracks
tracker=sim_tracks(tracker,tracker_par);
R=tracker_par.R;
V=zeros(1,size(z,2));
%first loop
for m=1:tracker.num_tracks
   % state propagation
   tracker.x(:,m) = tracker_par.F * tracker.x(:,m);
   P1 = tracker_par.F * tracker.P(:,:,m) * tracker_par.F' + tracker_par.Q;
   S = tracker_par.H*P1*tracker_par.H' + R;
   % zhat is the predicted measurement
   zhat=tracker_par.H*tracker.x(:,m);
   zi{m}=z-zhat*ones(1,size(z,2)); %innovation
   d=zi{m}'/S;
   zi_mah=abs(sum(zi{m}.*d',1));
   %measurement validation
   % measurement validation (using the validation gate, vg) 
   ng{m} = find(sqrt(zi_mah) <= tracker_par.vg);
   %mark shared measurements and keep number of times the meas are shared
   V(ng{m})=V(ng{m})+1; 
   pzz{m} = 1;
   %a-priori measurement probability
   s(ng{m})=tracker_par.pd*tracker_par.pw*tracker.PE(m)*pzz{m}/sum(pzz{m});
 end
 %second loop
 for n=1:length(V)
        v(n)=0;
        for k=1:V(n)
            v(n)=v(n)+s(n)./(1-s(n));
        end    
  end
  %third loop
  for m=1:tracker.num_tracks
        %probability of measurement visibility 
       pv=(1./(1-s(ng{m})))./(1+v(ng{m}));
       pv(isnan(pv))=0;
       %calculate probability of false alarm
       p=pv.*pzz{m}./tracker_par.pfa;
       delta_ipdaf = tracker_par.pd * tracker_par.pw * (1-sum(p));
       %propagate probabilities
       PM=tracker.PE(m);
       PM_hat=(1-delta_ipdaf)*PM/(1-delta_ipdaf*PM);
       beta0=(1-tracker_par.pd*tracker_par.pw)/(1-delta_ipdaf);    
       beta=tracker_par.pd*tracker_par.pw*p/(1-delta_ipdaf);
	   %combined inovation
	   zi_comb=[0];k_comb=zeros(1);
	   z1=z(:,ng{m});
	   zi1=zi{m}(:,ng{m});
	   for j=1:length(ng{m})
		   zi_comb=zi_comb+beta(j)*zi1(:,j);
		   k_comb=k_comb+beta(j)*zi1(:,j)*zi1(:,j)';
	   end
       %kalman gain
       K=P1*tracker_par.H'/S;
       %state estimate update
       tracker.x(:,m)=tracker.x(:,m)+K*zi_comb;
       %error covariance matrix update
       tracker.P(:,:,m)=P1-(1-beta0)*K*S*K'+K*(k_comb-zi_comb*zi_comb')*K';
       %markov chain update
       tracker.PE(m)=tracker_par.pmarkov*[PM_hat; 1]; 
end
warning on
