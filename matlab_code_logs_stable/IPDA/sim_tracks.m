function tracker=sim_tracks(tracker,tracker_par)
% SIM_TRACKS Deletion of similar tracks.
%    TRACKER = SIM_TRACKS(TRACKER, TRACK_THRES) Looks for similar tracks stored
%    in TRACKER and, based on the threshold value stored in the struct
%    TRACK_THRES, deletes all but one of the similar tracks.  For instance,
%    if there are two sets similar tracks: [1,3,4] and [2,5] then tracks
%    3, 4, and 5 would all be deleted.  The state vectors of tracks and various
%    other fields in TRACKER are updated and stored in the output argument.

a = find_sim_cols(tracker.x, tracker_par.sim_th);
tracker.x(:,a)=[];
tracker.P(:,:,a)=[];
tracker.PE(a)=[];

tracker.X_conf(:,a,:)=[];
tracker.X_tent(:,a,:)=[];
tracker.Cflag(a)=[];
%assign a brand new label to each additional track
L=tracker.new_lab;
N=size(tracker.x,2)-tracker.num_tracks;
tracker.lab(tracker.num_tracks+(1:N))=L+(1:N);
tracker.num_tracks = size(tracker.x, 2);
tracker.new_lab=L+N;