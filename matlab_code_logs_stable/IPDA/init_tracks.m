function tracker=init_tracks(z,z_old,tracker_par,tracker,tracker_init_val)
% INIT_TRACKS Initialization and augmentation of tracks.
%    TRACKER = INIT_TRACKS(z, Z_OLD,  TRACKER_PAR, TRACKER, ...
%                          TRACK_INIT_VAL) the current scene measurements, Z,
%    the previous scene measurements, Z_OLD, the structs containing the
%    clutter mapper, CLUTTERMAPPER_PAR, the tracker parameters, TRACKER_PAR, tracker
%    variable, TRACKER, and the tracker initialization values, TRACK_INIT_VAL.
%    The last argument, T, is an integer denoting the current time instant.
%    The augmented tracker is output to TRACKER.
%
%    For the definitions of the various structs above, see INIT_DATA,
%    INIT_PARS, INIT_VARS, INIT_SCENE.
%
% Sofia Suvorova
% Du Huynh, Oct 2002.
t=tracker.t;
[a,b] = find_sim_cols(z_old, z, tracker_par.vel_th);
a=a(:)';b=b(:)';
ydot=(z(b)-z_old(a))/tracker_par.delta;
y=[z(b);ydot];
siz = size(y,2);
num_tracks = tracker.num_tracks;
%probability of track existence 
tracker.PE(num_tracks+(1:siz))=tracker_par.pd*tracker_init_val.PE;


tracker.x = [tracker.x y];
for n=1:siz
   tracker.P(:,:,n+num_tracks) = tracker_init_val.P;
end
tracker.X_conf(:,num_tracks+(1:siz),t) = zeros(2,siz);
tracker.X_tent(:,num_tracks+(1:siz),t) = zeros(2,siz);
tracker.Cflag(num_tracks+(1:siz))=0;
