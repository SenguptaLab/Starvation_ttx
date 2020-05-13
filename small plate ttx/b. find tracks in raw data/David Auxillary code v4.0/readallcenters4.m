function tracks = readallcenters4(data_path,f1,f2)

% warning: all text files ni the folder must be centers!
%
% General: 
%   1. read all center files, 
%      but only rows corresponding to f1 < frame # < f2; 
%   2. track particles using track.m;  
%   3. arranges tracks data in tr structure format using putinfields.m
%     (no noise filtering apart from track.m) 
%
% Centers file format: 
% --------------------
%   1. TWO rows for every frame taken: 
%      X-coordinates of particles in frame followed by 
%      Y-coordinate of same particles in same frame
%   2. All coordinates are positive and in pixel units
%   3. Rows are filled with -9's for uniform length
%   4. First two rows of centers file are always all -9's
%   5. For some odd reason the LabView vi outputs the frame data 
%      in reverse order (FILO), so the last two rows are frame #1, 
%      second to last two rows are frame #2 etc. 
%
% Out format:
% -----------
% Two columns - [X Y] coordinates of particles 
% Third column - Frame number for each particle
% 
% tr = track(out, maxdisp, mem) input/output format:
%------------------------------------------------------
% out: as above.
% maxdisp: an estimate of the maximum distance (in pixels) that a particle 
%          would move in a single time interval.
% mem: this is the number of time steps that a particle can be
%      'lost' and then recovered again.  If the particle reappears
%      after this number of frames has elapsed, it will be
%      tracked as a new particle. The default setting is zero.
%      this is useful if particles occasionally 'drop out' of
%      the data.
% tr:  a list containing the original data rows sorted 
%      into a series of trajectories.  To the original input 
%      data structure there is appended an additional column 
%      containing a unique 'id number' for each identified 
%      particle trajectory.  The result array is sorted so 
%      rows with corresponding id numbers are in contiguous 
%      blocks, with the time variable a monotonically
%      increasing function inside each block.  For example:
%      For the input data structure:
%               (x)	        (y)	          (f)
%        out = 3.6000       5.0000      0.00000
%             15.1000      22.6000      0.00000
%              4.1000       5.5000      1.00000	
%    	      15.9000      20.7000      2.00000
%              6.2000       4.3000      2.00000
%  
% tr = track(pos,5,mem=2)
%  
% track will return the result 
%      		    (x)	        (y)	         (f) 	      (id)
%  		 tr = 3.60000      5.0000      0.00000      0.00000
%  		      4.10000      5.5000      1.00000      0.00000
%  		      6.20000      4.3000      2.00000      0.00000
%  		     15.1000      22.6000      0.00000      1.00000
%  		     15.9000      20.7000      2.00000      1.00000
%  
%  		for t=1 in the example above, one particle temporarily
%  		vanished.  As a result, the trajectory id=1 has one time
%  		missing, i.e. particle loss can cause time gaps to occur 
%       in the corresponding trajectory list. In contrast:
% tr = track(pos,5)
%  
% track will return the result 
%  			   (x)	        (y)	         (t) 	     (id)
%  		tr = 15.1000      22.6000      0.00000      0.00000
%             3.6000       5.0000      0.00000      1.00000
%             4.1000       5.5000      1.00000      1.00000
%             6.2000       4.3000      2.00000      1.00000
%            15.9000      20.7000      2.00000      2.00000
% 	
% 		where the reappeared 'particle' will be labelled as new
% 		rather than as a continuation of an old particle since
% 		mem=0.  It is up to the user to decide what setting of 
% 		'mem' will yeild the highest fidelity tracking.
%
% traque = putinfields(tr,1) and tracks (OUTPUT) format:
% ------------------------------------------------------
% an array of structures created by putinfields.m 
% Each structure corresponds to one track of one particle that was tracked
% Fields: tracks(i).x = x coordinates in pixels for track of particle i
%         tracks(i).y = y coordinates in pixels for track of particle i
%         tracks(i).f = corresponding frame numbers for track of particle i 
%         tracks(i).num = total number of frames in track of particle i
%         tracks(i).plate_id = number plate identifier for the plate of origin of particle i
%         otracks(i).min_edge = pixel value of edge with minimal Temperature
%         tracks(i).max_edge = pixel value of edge with maximal Temperature
%         tracks(i).min_T = Temperature at edge with minimal Temperature
%         tracks(i).max_T = Temperature at edge with maximal Temperature
%         tracks(i).total_time = total imaging time in seconds
%

global MAX_DIST 
global PARAM_MEM
global PARAM_DIM
global PARAM_GOOD
global PARAM_QUIET
global FPS
global SKIP_EVEN_LINES

PARAM.mem = PARAM_MEM;
PARAM.dim = PARAM_DIM;
PARAM.good = PARAM_GOOD;
PARAM.quiet = PARAM_QUIET;
tracks=[]; % 

warning off MATLAB:divideByZero;

files=dir(data_path);

plate_id = 0;
for i = 1:length(files) % go over all plates
    fn = files(i).name; 
    if length(fn)>11
        if strcmp(fn(end-10:end),'centers.txt')
            display(fn);
            bigmat = load([data_path filesep fn]);
            [rows,cols] = size(bigmat)

            if SKIP_EVEN_LINES % emulate 0.5 FPS data on 1 FPS file by skipping even line-pairs
                if nargin==3 && FPS==1 % first loop repetition ==> f1, f2, FPS need correction 
                    f1 = ceil(f1/2)
                    f2 = floor(f2/2)
                    FPS = 0.5
                end;
                tmp = [];
                for cur_row = rows:-4:2
                    tmp = [bigmat(cur_row-1,:); bigmat(cur_row,:); tmp];
                end;
                bigmat = tmp;
                [rows,cols] = size(bigmat)
            end;
            
            disp(['Frames: ' num2str(rows/2-1)]); % don't count two -9's rows
            out=[]; % three columns: X Y frame-number
            if nargin==1 % use default start/end times
                f1 = 1;
                f2 = rows/2-1;  % default: last row is all -9s
            elseif nargin ~= 3 % error: frame range not properly defined
                return;
            end;
            total_time = (f2-f1+1)/FPS;
            row_pair1 = rows/2-f1;
            row_pair2 = rows/2-f2;
            disp(['Actual time: ' num2str(total_time) ' (sec).']);
            for j=row_pair1:-1:row_pair2   % list out is in reverse order: 
                           % so frame 1 is last in bigmat.. 
                thisX = bigmat(2*j-1,:); % all x-coordinates at current frame (including -9s)
                thisY = bigmat(2*j,:);   % all y-coordinates at current frame (including -9s)
                thisF = rows/2-j; % REAL frame number
                choose=find(thisX > 0); % find particles (as opposed to -9's)
                if length(choose)>0 % add data from current frame to out 
                    out=[out;[thisX(choose)',thisY(choose)',thisF+zeros(length(choose),1)]];
                    % thisF+zeros(length(choose),1) is an array with
                    % length(choose) cells, all set to the value of the
                    % current real frame number. It will be 
                    % added as the third column for each particle
                end;
            end;
            tr = track(out, MAX_DIST, PARAM); % array of particle tracks for current centers file \
            tr = extrapolate_missing_frames4(tr);
            plate_id = plate_id+1; % individual plate identifier  
            %%% read cal file data %%%
            fn = regexprep(fn,'centers.txt','cal.csv');
            cal_data = load([data_path filesep fn]);
            fn_id = regexprep(fn,'cal.csv','');
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            traque = putinfields4(tr,1,plate_id,cal_data,total_time,fn_id,data_path);  
                                         % array of structures of individual particle tracks
                                         %            for current centers file
                                         % "1" value for clean flag -->
                                         % clean data from "bad" tracks
            tracks = [tracks, traque];  % add current structures of particle tracks to previous ones
        end; % if this is centers file
    end; % if length(fn) > 11
end; % i = 1:length(files) % go over all plates

    
    