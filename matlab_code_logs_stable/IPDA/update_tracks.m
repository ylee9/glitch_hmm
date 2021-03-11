 function tracker = update_tracks(tracker,tracker_par,t)
% UPDATE_TRACKS Update of tracks and optionally plotting of the scene and tracks.
%    TRACKER = UPDATE_TRACKS(TRACKER, TRACKER_PAR, T) Updates the confirmed,
%    tentative tracks stored in TRACKER and deletes tracks which have low
%    probabilities of detection.
%

% update confirmed, tentative tracks, etc
M=[];
PE=tracker.PE;
Cflag=tracker.Cflag;
for m=1:tracker.num_tracks
    tracker.X_conf(:,m,t)=[0;0];
    tracker.X_tent(:,m,t)=[0;0];
    x=tracker.x(:,m);
    if PE(m)<tracker_par.pd_low & Cflag(m)
        % kill previously existing tracks
        M=[M m];
    else
     if Cflag(m) | (~Cflag(m) & PE(m)>=tracker_par.pd_high)
        tracker.Cflag(m)=1;
        %confirmed track: any track with high probability of existence  
        tracker.X_conf(:,m,t)=x(1:2);
     else
        %tentative track: all and  only new tracks are tentative  
        tracker.X_tent(:,m,t)=x(1:2); 
     end
	end
end
   
% kill tracks with low prob of existance
if ~isempty(M)
    tracker.x(:,M)=[];
    tracker.P(:,:,M)=[];
    tracker.PE(M)=[];
    tracker.X_tent(:,M,:)=[];
    tracker.X_conf(:,M,:)=[];
    tracker.lab(M)=[];
    tracker.Cflag(M)=[];
    tracker.num_tracks=size(tracker.x, 2);
end
n=find(tracker.lab<=length(tracker.COUNT)); %existing tracks
m=find(tracker.lab>length(tracker.COUNT));  % new tracks
tracker.COUNT(tracker.lab(m))=1;
tracker.COUNT(tracker.lab(n))=tracker.COUNT(tracker.lab(n))+1;
%more cosmetics
if length(tracker.lab)~=tracker.num_tracks 
       tracker.lab=tracker.lab(1:tracker.num_tracks);
end
% true=sum(tracker.COUNT(tracker.lab)>=tracker_par.conf_length & tracker.Cflag);

    
