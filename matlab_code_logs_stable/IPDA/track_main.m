function [tracker]= track_main(meas,delta)
% TRACK_MAIN The main loop of target tracker.
% tracker parameters
H=[1 0];
sim_th=1e-5;
vel_th=1e-5;
pd_high=0.9;
pd_low=0.1;
R=1e1;
pfa=1e-6;
tracker_par = struct(...
	'delta',[],...
   'pmarkov', [0.96 0.02], ...               % Markov chain.
   'F', [], ...                               % State transition matrix.
   'H', H, ...
   'Q', [], ...                               % Covariance matrix of noise on the tracks.
   'R',R,...                                % measurement cov
   'pfa',pfa,...
   'vg', [], ...                            % Validation gate.
   'pd', 1, ...                            % Probability of detection
   'pw', 0.99, ...                           % Probability of being in the validation gate.
   'conf_length',0,...
   'pd_low', pd_low, ...                        % prob. of detection < pd_low => no tracks.
   'pd_high', pd_high, ...    % Prob. of detection > pd_high => confirmed tracks.% Between pd_low and pd_high => tentative tracks.
   'sim_th', sim_th, ...     % similarity threshold (Euclidean distance) between two tracks
   'vel_th', vel_th );  % velocity threshold for initialization of tracks
tracker_init_val = struct(...
   'P', eye(2)*0.0001, ...                    % initial value for each P cell array 
   'PE', 0.1,...        %PE 
    'q',1e-32 );   % fdot variance                          

tracker = struct(...
   't', 2, ...                  % time
   'P', [], ...                 % error covariance associated with the state
   'x', [], ...                 % mean the mth track
   'PE',[] , ...                % probability of tracks existence 
   'X_conf', [], ...            % confirmed tracks
   'X_tent', [], ...            % tentative tracks
   'new_lab',0,...              % label, which has to be assigned to a new track  
   'lab',[],...                 % track labels
   'Cflag',[],...               % confirmation flag (1 if track confirmed)   
   'COUNT',[],...               % counts how many times track label was displayed on screen
   'num_tracks',0 );           % current number of tracks

Nmeas=length(meas);
t = 2;
z=meas{1};
while t <= Nmeas
   z_old = z;
   z = meas{t};
   [tracker_par.F,tracker_par.Q] = comp_noise_cov(delta(t),tracker_init_val.q);
   tracker_par.delta=delta(t);
   tracker_par.vg=1/delta(t);
   tracker = LIPDAF(z_old,z,tracker_par,tracker,tracker_init_val);
   tracker = update_tracks(tracker,tracker_par,t);
   tracker.t = tracker.t + 1;
   t = tracker.t
end


