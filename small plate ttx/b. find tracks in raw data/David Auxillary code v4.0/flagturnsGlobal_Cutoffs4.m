function final_tracks=flagturnsGlobal_Cutoffs4(tr)
%
% General: 
%    1. Expand the tr structure to contain more information fields
%       (turns, durations, angles etc.)
%    2. Filter out what looks like noise tracks rather than worm tracks. 
%
% The goal here is to assemble the structure real tracks, to flag the
% points at which the worms turned and save info about the "runs" between 
% consecutive turns. 
% For each track of each worm the inifo includes the duration and 
% the length of individual "runs", the angle the worm was going at
% every point and on average during each tun, the indices at which the
% worms turned, etc.
% 
% INPUT: 
%  tr - an array of structures created by putinfields.m 
%       Each structure corresponds to one track of one worm 
%       Fields: tr(i).x = x coordinates in pixels for track of particle i
%               tr(i).y = y coordinates in pixels for track of particle i
%               tr(i).f = corresponding (real) frame numbers for track of particle i 
%               tr(i).num = total number of frames in track of particle i
%
% PARAMETERS: 
%  MIN_LEN - distance in pixels between points on track that are used to
%            determine velocity and angle (direction) of movements ; 
%            usually 1 pixel ~ 0.2 micron --> 6 pixels ~ worm body length 
%  MOVE_BY_FOR_TURN - minimal dx+dy total length in pixels to flag a turn
%  FLAG_TURNS3_N - number of points (N) used for calculating angles, directions etc.
%                  N-1 points at the start and N-1 points at the end of each
%                  track get ignored! N=5 was optimal according to Alex. 
%  MAX_FRAMES - Maximal number of frames used to determine angles
%  ANG_FOR_TURN - minimal turn angle to flag as "a turn"
%
% OUTPUT: 
%  tracks - an array of structures with newly added fields 
%           Each structure corresponds to one track of one worm
%           Old Fields from putinfields.m: 
%                   tracks(i).x = x coordinates in pixels at each frame for track of particle i
%                                 (column vector)
%                   tracks(i).y = y coordinates in pixels at each frame for track of particle i
%                                 (column vector)
%                   tracks(i).f = corresponding frame numbers for track of particle i 
%                                 (column vector)
%                   tracks(i).num = total number of frames in track of particle i
%                   tracks(i).plate_id = plate id for particle i
%            Newly added:
%                   tracks(i).t = corresponding (real) times (in seconds) for track of particle i 
%                 Segments: 
%                   tracks(i).seg_indx = segs;
%                   tracks(i).seg_ang = seg_ang;
%                   tracks(i).seg_num = Ns; 
%                   tracks(i).seg_dx = seg_dx; 
%                   tracks(i).seg_dy = seg_dy;
%                   tracks(i).seg_dist = sqrt(seg_dx.^2+seg_dy.^2);
%                   tracks(i).seg_dur = seg_dt;
%                 Runs:
%                   tracks(i).run_duration = durations of individual runs in number of frames
%                   tracks(i).run_time = time in seconds of individual runs
%                   tracks(i).run_dx = end to end X-distance in pixels for every run
%                   tracks(i).run_dy = end to end Y-distance in pixels for every run
%                   tracks(i).run_dist = end-to-end distance (NOT along curve) in pixels for every run
%                   tracks(i).run_Int_along_path = distance along curve, in pixels, for each run 
%                   tracks(i).curved = 1 if run_Int_along_path > pi/2*run_dist; 0 otherwise
%                   tracks(i).run_angle = average directions of individual runs in radians
%                                         -pi < angle <= pi , 
%                                         pi/2=up, 0 p/m epsilon=right, -pi/2=down, pi or -pi=left
%                   tracks(i).run_indx = index pairs ([s1,e1; s2,e2; ..]) of run starts and run
%                                        ends in vectors such as x,y,t.
%                                        E.G., t1=t(run_indx(Nr,1)) = starting
%                                              time of the Nr-th run
%                                              t2=t(run_indx(Nr,2)) = ending
%                                              time of the Nr-th run
%                                              (t2-t1) = duration of Nr-th run
%                   tracks(i).run_num = total number of runs
%                   tracks(i).turn_indx = indices of turns between runs
%
%
%               Isothermal Tracks:
%                   tracks(i).IT_duration = durations of individual Isothermal Tracks in number of frames
%                   tracks(i).IT_time = time in seconds of individual ITs
%                   tracks(i).IT_dx = end to end X-distance in pixels for every IT
%                   tracks(i).IT_dy = end to end Y-distance in pixels for every IT
%                   tracks(i).IT_dist = total distance (along curve) in pixels for every IT
%                   tracks(i).IT_end2end_dist = distance in pixels in straight line from start to end for every IT
%                   tracks(i).IT_angle = average directions of individual ITs in radians
%                                         -pi < angle <= pi , 
%                                         pi/2=up, 0 p/m epsilon=right, -pi/2=down, pi or -pi=left
%                   tracks(i).IT_indx = index pairs ([s1,e1; s2,e2; ..]) of IT starts and IT
%                                        ends in vectors such as x,y,t.
%                                        E.G., t1=t(run_indx(n,1)) = starting
%                                              time of the n-th IT
%                                              t2=t(run_indx(n,2)) = ending
%                                              time of the n-th IT
%                                              (t2-t1) = duration of n-th IT
%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
global FPS

global SEG_D      % distance in pixels that counts as a discrete "step"
global SEG_T1     % minimal time in seconds for a segment ro be counted
global SEG_T2     % maximal time in seconds for a segment ro be counted
global N_SEG      % minimal number of segments for "run"  
global ALPHA      % minimal changle in angle that counts for a turn
global MAXDANG_IT % max deviation from pi/2 (or -pi/2) which is still considered IT (radians)
global MINTIME_IT % min time of track to be considered IT (in seconds)

sim_ang = inline('abs(a-b)<c | ((pi-abs(a)<c/2) & (pi-abs(b)<c/2))','a','b','c');
% -pi < each angle <= pi , 
% pi/2=up, 0 p/m epsilon=right, -pi/2=down, pi or -pi=left

tracks=tr;
for i=1:length(tr) %for each track of each worm
    x = tracks(i).x;
    y = tracks(i).y;
    f = tracks(i).f;
    N = tracks(i).num; % length of x,y,f
    
    t = f / FPS; 
    tracks(i).t = t;

    %%%%%%%%%%%%% plot times of tracks %%%%%%%%%%%%%
    %figure(1000+tracks(i).plate_id);
    %hold on;
    %plot([t(1), t(end)],[i,i]);
    %title(['Tracks found for plate no. ' num2str(tracks(i).plate_id)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    dist=[0; cumsum(sqrt(diff(x).^2+diff(y).^2))]; % cummulative distance in pixels
    
    tmp_segs = [1]; % indices delimiting potential segments
    for k=2:N       % (not taking into account time criteria yet)
        if (dist(k) > SEG_D)
            dist = dist-dist(k); 
            tmp_segs = [tmp_segs; k];
        end;
    end;
    if tmp_segs(end)<N
        tmp_segs = [tmp_segs; N];
    end;
    
    segs = [];    % pairs of indices delimiting real segments 
    seg_ang = []; % angle of each segment
    seg_dx = [];
    seg_dy = [];
    seg_dt = [];
    % e.g., 
    % segs     = [[1,10]; [11,20]; [23,31]; [35,44]; [46,53]; [55,62]; ...]
    % seg_ang  = [ 0.1;    0.12;    0.09;    1.05;     1.12;   1.97; ...]
    
    for k = 1:(length(tmp_segs)-1)
        dt = t(tmp_segs(k+1))-t(tmp_segs(k));
        if (dt>SEG_T1) & (dt<SEG_T2) % time criteria for real segments
            segs = [segs; [tmp_segs(k), tmp_segs(k+1)]];
        end;
    end;
    Ns = size(segs,1);
    
    for k = 1:Ns
        dx = x(segs(k,2))-x(segs(k,1));
        dy = y(segs(k,2))-y(segs(k,1));
        dt = t(segs(k,2))-t(segs(k,1));
        seg_ang = [seg_ang; angle(dx+dy*sqrt(-1))];
        seg_dx = [seg_dx; dx];
        seg_dy = [seg_dy; dy];
        seg_dt = [seg_dt; dt];
        % -pi < seg_ang(:) <= pi , 
        % pi/2=up, 0 p/m epsilon=right, -pi/2=down, pi or -pi=left
    end;
    
    tracks(i).seg_indx = segs;
    tracks(i).seg_ang = seg_ang;
    tracks(i).seg_num = Ns; 
    tracks(i).seg_dx = seg_dx; 
    tracks(i).seg_dy = seg_dy;
    tracks(i).seg_dist = sqrt(seg_dx.^2+seg_dy.^2);
    tracks(i).seg_dur = seg_dt;
    % didn't save seg frames or real time - but it is easy to do here
    
    % Group segments into groups that will become "runs"
    % e.g., 
    % segs       = [[1,10]; [11,20]; [23,31]; [35,44]; [46,53]; [55,62]; ...]
    % seg_ang    = [ 0.1;    0.12;    0.09;    1.05;     1.12;    1.97; ...]
    % run_indx   = [[1,31]; [35,53]; [55,..] ...]
    % run_ang    = [ 0.11;   1.09;   ...]
    % turns_indx = [  33;     54;    ...]
    gr_id = 1;     % group identifier
    cur_gr = [1];  % current group
    segs(1,3) = gr_id;
    for k = 2:Ns % group segments in "similar direction" groups
        cur_gr_ang = mean(seg_ang(cur_gr));  % current group angle (mean of group)
        
        %%% Continue runs as long as there are no changes to AVERAGE direction %%%
        continue_run = sim_ang(cur_gr_ang, seg_ang(k), ALPHA); % two angles not closer than ALPHA
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if continue_run
            cur_gr = [cur_gr, k]; % add k-th segment to current group
                                  % do not change group id
        else % new segment angle is different than group by more than ALPHA 
            segs(cur_gr,4) = cur_gr_ang; % save mean group angle in 4th column
            cur_gr = [k];            % make new group, starting with segment k
            gr_id = gr_id+1;         % make new group id
                                     % (new group angle will be calculated
                                     %  in next loop iteration) 
        end;
        segs(k,3) = gr_id; % add group identifier to 3rd column of segs        
    end;
    
    run_indx = [];
    run_ang = [];
    run_dur = [];
    run_dx = [];
    run_dy = [];
    run_dist = [];
    run_Int_along_path = [];
    turn_indx = [];
    run_time = [];
    Nr = 0;
    cur_gr = 1;
    cur_first = 1;
    cur_gr_size = 0;
    for k = 1:Ns % find groups that are larger than N_SEG runs, calculate angles.. 
        if segs(k,3)==cur_gr
            cur_gr_size = cur_gr_size+1;
        else     % save current run, start new run
            if cur_gr_size >= N_SEG % "run" includes segments (cur_first)..(k-1), all in cur_gr
                Nr = Nr+1;         % segment (k) may be the first segment in the next run
                run_indx(Nr,1) = segs(cur_first,1); % index of start of segment (cur_first)
                run_indx(Nr,2) = segs(k-1,2);       % index of end   of segment (k-1)
                dx = x(run_indx(Nr,2))-x(run_indx(Nr,1));
                dy = y(run_indx(Nr,2))-y(run_indx(Nr,1)); 
                df = f(run_indx(Nr,2))-f(run_indx(Nr,1));
                dt = t(run_indx(Nr,2))-t(run_indx(Nr,1));
                run_dx  = [run_dx;  dx]; % pixels
                run_dy  = [run_dy;  dy]; % pixels
                run_dur = [run_dur; df]; % frames
                run_time = [run_time; dt]; % seconds
                run_dist = [run_dist; sqrt(dx^2+dy^2)]; % end to end distance
                                                        % (not intergal along path)
                tmp = diff(x(run_indx(Nr,1):run_indx(Nr,2))).^2; 
                tmp = tmp+diff(y(run_indx(Nr,1):run_indx(Nr,2))).^2; 
                run_Int_along_path = [run_Int_along_path; sum(sqrt(tmp))]; % intergal along path
                run_ang = [run_ang; angle(dx+dy*sqrt(-1))];
                % -pi < seg_ang(:) <= pi , 
                % pi/2=up, 0 p/m epsilon=right, -pi/2=down, pi or -pi=left
            end;
            cur_gr = segs(k,3);
            cur_first = k;
            cur_gr_size = 1; 
        end;  % else: save current run
    end; % for k=... k has reached Ns
    if cur_gr_size >= N_SEG % add last run (wasn't counted because loop ended)
       Nr = Nr+1;
       run_indx(Nr,1) = segs(cur_first,1); % index of start of segment (cur_first)
       run_indx(Nr,2) = segs(k-1,2);       % index of end   of segment (k-1)
       dx = x(run_indx(Nr,2))-x(run_indx(Nr,1));
       dy = y(run_indx(Nr,2))-y(run_indx(Nr,1)); 
       df = f(run_indx(Nr,2))-f(run_indx(Nr,1));
       dt = t(run_indx(Nr,2))-t(run_indx(Nr,1));
       run_dx  = [run_dx;  dx]; % pixels
       run_dy  = [run_dy;  dy]; % pixels
       run_dur = [run_dur; df]; % frames
       run_time = [run_time; dt]; % seconds
       run_dist = [run_dist; sqrt(dx^2+dy^2)]; % end to end distance
                                               % (not intergal along path)
       tmp = diff(x(run_indx(Nr,1):run_indx(Nr,2))).^2; 
       tmp = tmp+diff(y(run_indx(Nr,1):run_indx(Nr,2))).^2; 
       run_Int_along_path = [run_Int_along_path; sum(sqrt(tmp))]; % intergal along path 
       run_ang = [run_ang; angle(dx+dy*sqrt(-1))];
       % -pi < seg_ang(:) <= pi , 
       % pi/2=up, 0 p/m epsilon=right, -pi/2=down, pi or -pi=left
    end; % add last run
 
    tracks(i).run_indx = run_indx; % double column vector of start/end indices for each run
    tracks(i).run_duration = run_dur; % column vector of durations (in seconds) for each run
    tracks(i).run_angle = run_ang; % column vector of angles for each run
    tracks(i).run_dx = run_dx; % column vector of total x-axis distance (in pixels) for each run
    tracks(i).run_dy = run_dy; % column vector of total y-axis distance (in pixels) for each ru
    tracks(i).run_dist = run_dist;  % column vector of end-to-end distance (in pixels) for each ru
    tracks(i).run_Int_along_path = run_Int_along_path; % column vector of integrated-along-curve distance (in pixels) for each ru
    tracks(i).curved = run_Int_along_path > (pi/2 * run_dist); % column vector, Alex/Damon criteria for curved run
    tracks(i).run_num = Nr;
    tracks(i).run_time = run_time;
    %%% save Isothermal Tracks %%%
    IT_index = find( (abs(abs(run_ang)-pi/2) < MAXDANG_IT) & (run_time > MINTIME_IT) );
    tracks(i).IT_indx = run_indx(IT_index,:); 
    tracks(i).IT_duration = run_dur(IT_index);
    tracks(i).IT_angle = run_ang(IT_index);
    tracks(i).IT_dx = run_dx(IT_index);
    tracks(i).IT_dy = run_dy(IT_index);
    tracks(i).IT_dist = run_dist(IT_index);
    tracks(i).IT_time = run_time(IT_index);
    tracks(i).IT_num = length(IT_index);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for k=1:(Nr-1) % calculate index vector of turns (between runs)
        indx1 = run_indx(k,2); % end of k-th run
        indx2 = run_indx(k+1,1); % start of (k+1)-th run
        turn_indx = [turn_indx; round(mean([indx1, indx2]))];
    end;
    tracks(i).turn_indx = turn_indx;  
    
    %display([num2str(Nr) ' runs,   ' num2str(length(turn_indx)) ' turns']);
    
end; % for i=1:length(tr) , i.e., for each track of each worm

%%% Impose Alex-Damon cut-offs %%%
%%% tracks must be longer than one minute (total) %%%
%%% runs that start at the first frame of a track are ignored %%%
%%% runs that end on the last frame of a track are ignored %%%
k = 0;
for i=1:length(tracks)  
    if (max(tracks(i).t)-min(tracks(i).t))>60  % longer than 1 min
        tmp = tracks(i);
        keep_this = 1;
        if size(tmp.run_indx,1)<1 || size(tmp.run_indx,2)<1
            keep_this = 0;
        end;
        if keep_this && tmp.run_indx(1,1)==1 % get rid of first run 
            if tmp.run_num < 2 
                keep_this = 0; 
            else
                tmp.run_indx = tmp.run_indx(2:end,:);
                tmp.run_dx  = tmp.run_dx(2:end,:);
                tmp.run_dy  = tmp.run_dy(2:end,:);
                tmp.run_duration = tmp.run_duration(2:end,:);
                tmp.run_time = tmp.run_time(2:end,:);
                tmp.run_dist = tmp.run_dist(2:end,:);
                tmp.run_Int_along_path = tmp.run_Int_along_path(2:end,:);
                tmp.run_angle = tmp.run_angle(2:end,:);
                tmp.run_num = tmp.run_num-1;
                tmp.curved = tmp.curved(2:end,:); 
            end;
        end; % get rid of first run
        if keep_this && tmp.run_indx(end,2)==length(tmp.f) % get rid of last run 
            if tmp.run_num < 2 
                keep_this = 0; 
            else
                tmp.run_indx = tmp.run_indx(1:(end-1),:);
                tmp.run_dx  = tmp.run_dx(1:(end-1));
                tmp.run_dy  = tmp.run_dy(1:(end-1));
                tmp.run_duration = tmp.run_duration(1:(end-1));
                tmp.run_time = tmp.run_time(1:(end-1));
                tmp.run_dist = tmp.run_dist(1:(end-1));
                tmp.run_Int_along_path = tmp.run_Int_along_path(1:(end-1));
                tmp.run_angle = tmp.run_angle(1:(end-1));
                tmp.run_num = tmp.run_num-1;
                tmp.curved = tmp.curved(1:(end-1)); 
            end;
        end; % get rid of last run 
        if keep_this
            k = k+1; 
            tracks2(k) = tmp;
        end;
    end;
end;
tracks = tracks2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

final_tracks = tracks; % filter_out_some_noise moved to analyze_runs_and_turns.m

return;
