function [data_summary]=make_data_summary4(tr)

global ANG_FOR_RW
NUM_INTERVALS_FOR_ARS=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set default values + comments describing the fields.           %%%
%%% Data structure:                                                %%%
%%% data_summary.plates.plate_id - cell array with prefix of       %%%
%%%                                plate file-names (string).      %%%
%%% data_summary.num_of_plates - total number of plates.           %%%
%%% Cell arrays - each element is a column vector with values      %%%
%%%               corresponding to one plate (e.g., run-times,     %%%
%%%               IT_temperatures, etc.)                           %%%
%%% Example of syntax: C = {[1;2;3],[4;5],[6;7;8;9]}               %%%
%%%               C{2} == [4,5] i.e., equivalent to cell2mat(C(2)) %%%
%%% Arrays - column vectors where each element is a numerical      %%%
%%%          value corresponding to one plate (e.g., Cryophilic    %%%
%%%          index, total number of runs down the gradient, etc.)  %%%
%%%                                                                %%%
%%% See detailed descriptions in comments below.                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_summary.path = [];               % path to directory with raw data
data_summary.mean_relative_error = 0; % mean value (over all runs) of the ratio 
                                      % (distance of raw data points from run-line)/(run-length)
% All runs from all plates 
data_summary.runs.time = {};     % duration of each run in sec (by plates)
data_summary.runs.dx = {};
data_summary.runs.angle = {};    % angle of each run in rad (by plates) 
data_summary.runs.velocity = {}; % average velocity during each run in pxl/sec (by plates)
data_summary.runs.velocity_mm = {};

% All Isothermal Tracks
data_summary.ITs.time = {};        % duration of each IT in sec (by plates)
data_summary.ITs.temperature = {}; % average temperature of each IT in deg C (by plates)
data_summary.ITs.velocity = {};    % average velocity during each IT in pxl/sec (by plates)
data_summary.ITs.velocity_mm = {};
data_summary.ITs.X_pos = {};
data_summary.ITs.start_end = {};
% All turns from all plates
data_summary.turns.t = {};          % timing in sec of each turn (by plates)
data_summary.turns.x = {};          % x pos in px of each turn (by plates)
data_summary.turns.y = {};          % y pos in px of each turn (by plates)
data_summary.turns.from_angle = {}; % angle of run preceding turn in rad (by plates)
data_summary.turns.to_angle = {};   % angle of run following turn in rad (by plates) 

% Data per plate 
data_summary.plates.plate_id = {};     % string identifying plate file-name
data_summary.plates.num_of_plates = 0; % total number of plates,  
data_summary.plates.av_velocity = [];  % average velocity per plate
data_summary.plates.av_velocity_mm = [];
data_summary.plates.max_T = [];
data_summary.plates.min_T = [];
data_summary.plates.max_edge = [];
data_summary.plates.min_edge = [];
data_summary.plates.reversed = [];


data_summary.plates.runs_total_time_dn = []; % sum durations of all runs down in sec (per plate)
data_summary.plates.runs_total_time_vr = []; % sum durations of all runs vertical in sec (per plate)
data_summary.plates.runs_total_time_up = []; % sum durations of all runs up in sec (per plate)
data_summary.plates.runs_total_time_all = []; % sum durations of all runs  in sec (per plate)
data_summary.plates.runs_total_time_hr = []; % sum durations of all runs up and down in sec (per plate)
data_summary.plates.runs_total_dx_dn = []; % sum dx of all runs down in sec (per plate)
data_summary.plates.runs_total_dx_vr = []; % sum dx of all runs vertical in sec (per plate)
data_summary.plates.runs_total_dx_up = []; % sum dx of all runs up in sec (per plate)
data_summary.plates.runs_N_dn = []; % number of all runs up in sec (per plate)
data_summary.plates.runs_N_vr = []; % number of all runs up in sec (per plate)
data_summary.plates.runs_N_up = []; % number of all runs up in sec (per plate)
data_summary.plates.runs_N_all = []; % number of all runs  in sec (per plate)
data_summary.plates.runs_N_hr = []; % number of all runs up and down in sec (per plate)

data_summary.plates.I_cr = []; % Cryophilic Index (per plate)
data_summary.plates.I_cr_dx = [];

data_summary.plates.ITs_total_time = []; % sum durations of all ITs in sec(per plate)
data_summary.plates.ITs_total_dist = []; % sum durations of all ITs in sec(per plate)
data_summary.plates.ITs_total_dist_mm = [];
data_summary.plates.ITs_N = [];          % number of all ITs (per plate)
data_summary.plates.ITs_frac_time = []; % fraction of time spent tracking(per plate)
data_summary.plates.ITs_frac_dist= []; % fraction of length spent tracking(per plate)
data_summary.plates.ITs_frac_N = []; % fraction of time spent tracking(per plate)
data_summary.plates.ITs_av_temp = [];    % average temperature of all ITs in deg C (per plate)
data_summary.plates.ITs_av_temp_pos = [];
data_summary.plates.ITs_temp_band_width = []; % with of IT band in deg C (per plate)
data_summary.plates.ITs_avg_length = []; % Avg length of IT in px
data_summary.plates.ITs_avg_length_mm = [];

data_summary.plates.px2mm = []; % approx conversion from pixels to mm (mm/px)
data_summary.plates.turns_num_in_intrvl = []; % number of turns in 5 min bins intervals (per plate)
                                              % each column vector corresponds to the number of   
                                              % turns in 1-5, 6-10, 11-15.. min for a given plate 
data_summary.plates.total_turns = []; % number of turns overall

% All plates combined  
% [mean-over-plates, standard-error-over-plates]
data_summary.ALL.all_mean_velocity = [];     % (pxl/sec)
data_summary.ALL.all_mean_velocity_mm = [];
data_summary.ALL.run_mean_velocity = [];     % (pxl/sec)
data_summary.ALL.run_mean_velocity_mm = [];
data_summary.ALL.run_dn_total_duration = []; % (sec)
data_summary.ALL.run_vr_total_duration = []; % (sec)
data_summary.ALL.run_up_total_duration = []; % (sec)
data_summary.ALL.run_all_total_duration = [];
data_summary.ALL.run_hr_total_duration = [];
data_summary.ALL.run_dn_mean_duration =[];
data_summary.ALL.run_vr_mean_duration =[];
data_summary.ALL.run_up_mean_duration =[];
data_summary.ALL.run_all_mean_duration =[];
data_summary.ALL.run_hr_mean_duration=[];
data_summary.ALL.run_dn_num = [];            % ()
data_summary.ALL.run_vr_num = [];            % ()
data_summary.ALL.run_up_num = [];            % ()
data_summary.ALL.run_all_num =[];
data_summary.ALL.run_hr_num =[];
data_summary.ALL.I_cr = [];                  % ()
data_summary.ALL.I_cr_dx = [];
data_summary.ALL.IT_mean_velocity = [];      % (pxl/sec)
data_summary.ALL.IT_mean_velocity_mm = [];
data_summary.ALL.IT_mean_duration = [];      % (sec)
data_summary.ALL.IT_mean_temperature = [];   % (deg C)
data_summary.ALL.IT_band_width = [];         % (deg C)
data_summary.ALL.IT_total_duration = [];     % (sec)
data_summary.ALL.IT_num = [];                % ()
data_summary.ALL.ITs_frac_time = [];          % (s/s n)
data_summary.ALL.ITs_frac_dist = [];          % ()
data_summary.ALL.ITs_frac_N = [];             % ()
data_summary.ALL.ITs_avg_length = [];
data_summary.ALL.ITs_avg_length_mm = [];
data_summary.ALL.turns_num_in_intrvl = [];   % 1st column: mean , 2nd column: standard error
data_summary.ALL.total_turns =[];            % ()   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MB added Units
data_summary.plates.units.plate_id = '';
data_summary.plates.units.num_of_plates = '';
data_summary.plates.units.av_velocity = 'px/s';
data_summary.plates.units.av_velocity_mm = 'mm/s';
data_summary.plates.units.runs_total_time_dn = 's';
data_summary.plates.units.runs_total_time_vr = 's';
data_summary.plates.units.runs_total_time_up = 's';
data_summary.plates.units.runs_total_time_all = 's';
data_summary.plates.units.runs_total_time_hr = 's';
data_summary.plates.units.runs_total_dx_dn = 'px';
data_summary.plates.units.runs_total_dx_vr = 'px';
data_summary.plates.units.runs_total_dx_up = 'px';
data_summary.plates.units.runs_N_dn = '';
data_summary.plates.units.runs_N_vr = '';
data_summary.plates.units.runs_N_up = '';
data_summary.plates.units.runs_N_all = '';
data_summary.plates.units.runs_N_hr = '';
data_summary.plates.units.I_cr = '';
data_summary.plates.units.I_cr_dx = '';
data_summary.plates.units.ITs_total_time = 's';
data_summary.plates.units.ITs_total_dist = 'px';
data_summary.plates.units.ITs_total_dist_mm = 'mm';
data_summary.plates.units.ITs_frac_time = 's/s no units';
data_summary.plates.units.ITs_frac_dist = 'px/px no units';
data_summary.plates.units.ITs_frac_N = 'no units';
data_summary.plates.units.ITs_N = '';
data_summary.plates.units.ITs_av_temp = 'C';
data_summary.plates.units.ITs_av_temp_pos = 'C';
data_summary.plates.units.ITs_temp_band_width = 'C';
data_summary.plates.units.ITs_avg_length = 'px';
data_summary.plates.units.ITs_avg_length_mm = 'mm';
data_summary.plates.units.px2mm = 'mm/px';
data_summary.plates.units.turns_num_in_intrvl = '';
data_summary.plates.units.total_turns = '';


data_summary.ALL.units.all_mean_velocity = 'px/s';
data_summary.ALL.units.all_mean_velocity_mm = 'mm/s';
data_summary.ALL.units.run_mean_velocity = 'px/s';
data_summary.ALL.units.run_mean_velocity_mm = 'mm/s';
data_summary.ALL.units.run_dn_total_duration = 's';
data_summary.ALL.units.run_vr_total_duration = 's';
data_summary.ALL.units.run_up_total_duration = 's';
data_summary.ALL.units.run_all_total_duration = 's';
data_summary.ALL.units.run_hr_total_duration = 's';
data_summary.ALL.units.run_dn_num = '';
data_summary.ALL.units.run_vr_num = '';
data_summary.ALL.units.run_up_num = '';
data_summary.ALL.units.run_all_num = '';
data_summary.ALL.units.run_hr_num = '';
data_summary.ALL.units.I_cr = '';
data_summary.ALL.units.I_cr_dx = '';
data_summary.ALL.units.IT_mean_velocity = 'px/s';
data_summary.ALL.units.IT_mean_velocity_mm = 'mm/s';
data_summary.ALL.units.IT_mean_duration = 's';
data_summary.ALL.units.IT_mean_temperature = 'C';
data_summary.ALL.units.IT_band_width = 'C';
data_summary.ALL.units.IT_total_duration = 's';
data_summary.ALL.units.IT_num = '';
data_summary.ALL.units.ITs_frac_time = 's/s no units';
data_summary.ALL.units.ITs_frac_dist = 'px/px no units';
data_summary.ALL.units.ITs_frac_N = 'no units';
data_summary.ALL.units.turns_num_in_intrvl = '';
data_summary.ALL.units.total_turns = '';
data_summary.ALL.units.run_dn_mean_duration = 's';
data_summary.ALL.units.run_vr_mean_duration = 's';
data_summary.ALL.units.run_up_mean_duration = 's';
data_summary.ALL.units.run_all_mean_duration = 's';
data_summary.ALL.units.run_hr_mean_duration = 's';
data_summary.ALL.units.ITs_avg_length = 'px';
data_summary.ALL.units.ITs_avg_length_mm = 'mm';
data_summary.ALL.units.px2mm = 'mm/px';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Real values from tr: %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_summary.path = tr(1).path; % path to directory with raw data
[p fl] = fileparts(tr(1).path);
data_summary.folder = fl;

run_err = [];
run_dist = [];
for k = 1:length(tr)
    run_err = [run_err; tr(k).run_err];
    run_dist = [run_dist; tr(k).run_dist];
end;
indx = find(run_dist>0); % should be all but just in case
data_summary.mean_relative_error = mean(run_err(indx)./run_dist(indx));


current_plate = []; 
cur_plt_time = 0;
cur_plt_distance = 0;
cur_plt_distance_mm = 0;
for k = 1:(length(tr)+1) % +1 so that last plate can be recorded conveniently 
    
    % record previous plate if necessary: 
    %------------------------------------
    rec_last_plate = (k==(length(tr)+1)); % last value of k --> only record last plate
    if ~rec_last_plate % check if finished with current (not last) plate
        rec_new_plate = ~isequal(current_plate, tr(k).file_name_id) && k>1; 
    end;
    
    if rec_new_plate || rec_last_plate % need to record a plate, last or not last
        display(['Recording plate: ' current_plate]);
        prev_plate_good = (cur_plt_time>0) && (~isempty(this_plate.runs.time)); % check if plate has data
        prev_plate_good = prev_plate_good &&  (~isempty(this_plate.turns.t));   % check if plate has data
        if  prev_plate_good % record data for previous plate
            % Runs
            data_summary.runs.time = [data_summary.runs.time, {this_plate.runs.time}];
            data_summary.runs.dx = [data_summary.runs.dx, {this_plate.runs.dx}];
        	data_summary.runs.angle = [data_summary.runs.angle, {this_plate.runs.angle}];
        	data_summary.runs.velocity = [data_summary.runs.velocity, {this_plate.runs.velocity}];
            data_summary.runs.velocity_mm = [data_summary.runs.velocity_mm, {this_plate.runs.velocity_mm}];
            % ITs    
         	data_summary.ITs.time = [data_summary.ITs.time, {this_plate.ITs.time}];       
         	data_summary.ITs.temperature = [data_summary.ITs.temperature, {this_plate.ITs.temperature}]; 
        	data_summary.ITs.velocity = [data_summary.ITs.velocity, {this_plate.ITs.velocity}];  
            data_summary.ITs.velocity_mm = [data_summary.ITs.velocity_mm, {this_plate.ITs.velocity_mm}];
            data_summary.ITs.X_pos = [data_summary.ITs.X_pos, {this_plate.ITs.X_pos}];
            data_summary.ITs.start_end =  [data_summary.ITs.start_end, {this_plate.ITs.start_end}];
            % Turns
         	data_summary.turns.t = [data_summary.turns.t, {this_plate.turns.t}];   
            data_summary.turns.x = [data_summary.turns.x, {this_plate.turns.x}];   
            data_summary.turns.y = [data_summary.turns.y, {this_plate.turns.y}];   
         	data_summary.turns.from_angle = [data_summary.turns.from_angle, {this_plate.turns.from_angle}]; 
         	data_summary.turns.to_angle = [data_summary.turns.to_angle, {this_plate.turns.to_angle}];                 
            % Plates (general, runs, ITs, turns)     
            data_summary.plates.plate_id = [data_summary.plates.plate_id, {current_plate}];
            data_summary.plates.av_velocity = [data_summary.plates.av_velocity; cur_plt_distance/cur_plt_time];
            data_summary.plates.av_velocity_mm = [data_summary.plates.av_velocity_mm; cur_plt_distance_mm/cur_plt_time];
            data_summary.plates.num_of_plates = data_summary.plates.num_of_plates + 1;  
            data_summary.plates.max_T = [data_summary.plates.max_T, maxT  ];
            data_summary.plates.min_T = [data_summary.plates.min_T, minT  ];
            data_summary.plates.max_edge = [data_summary.plates.max_edge, maxedge  ];
            data_summary.plates.min_edge = [data_summary.plates.min_edge, minedge  ];
            data_summary.plates.reversed = [data_summary.plates.reversed, reversed  ];

            
        	data_summary.plates.runs_total_time_dn = [data_summary.plates.runs_total_time_dn; this_plate.plates.runs_total_time_dn];
        	data_summary.plates.runs_total_time_vr = [data_summary.plates.runs_total_time_vr; this_plate.plates.runs_total_time_vr];
        	data_summary.plates.runs_total_time_up = [data_summary.plates.runs_total_time_up; this_plate.plates.runs_total_time_up];
        	data_summary.plates.runs_total_time_all = [data_summary.plates.runs_total_time_all; this_plate.plates.runs_total_time_all];
        	data_summary.plates.runs_total_time_hr = [data_summary.plates.runs_total_time_hr; this_plate.plates.runs_total_time_hr];
            
            data_summary.plates.runs_total_dx_dn = [data_summary.plates.runs_total_dx_dn; this_plate.plates.runs_total_dx_dn];
            data_summary.plates.runs_total_dx_vr = [data_summary.plates.runs_total_dx_vr; this_plate.plates.runs_total_dx_vr];
            data_summary.plates.runs_total_dx_up = [data_summary.plates.runs_total_dx_up; this_plate.plates.runs_total_dx_up];

            data_summary.plates.runs_N_dn = [data_summary.plates.runs_N_dn; this_plate.plates.runs_N_dn];
        	data_summary.plates.runs_N_vr = [data_summary.plates.runs_N_vr; this_plate.plates.runs_N_vr];
        	data_summary.plates.runs_N_up = [data_summary.plates.runs_N_up; this_plate.plates.runs_N_up];
            data_summary.plates.runs_N_all = [data_summary.plates.runs_N_all; this_plate.plates.runs_N_all];
        	data_summary.plates.runs_N_hr = [data_summary.plates.runs_N_hr; this_plate.plates.runs_N_hr];

        	I_cr = (this_plate.plates.runs_total_time_dn - this_plate.plates.runs_total_time_up);
            I_cr = I_cr/(this_plate.plates.runs_total_time_dn + this_plate.plates.runs_total_time_up);
            I_cr_dx = -(this_plate.plates.runs_total_dx_dn + this_plate.plates.runs_total_dx_up);
            I_cr_dx = I_cr_dx/(abs(this_plate.plates.runs_total_dx_dn) + abs(this_plate.plates.runs_total_dx_up));
            if(reversed)
                I_cr = -I_cr;
                I_cr_dx = -I_cr_dx;
            end
            
            
        	data_summary.plates.I_cr = [data_summary.plates.I_cr; I_cr];
            data_summary.plates.I_cr_dx = [data_summary.plates.I_cr_dx; I_cr_dx];
            

            if length(this_plate.ITs.time) > 0 % record mean IT temperature and Tracking Band width for previous plate
                Pr = this_plate.ITs.time ./ sum(this_plate.ITs.time);
                avT = sum(this_plate.ITs.temperature .* Pr);
                BW = sqrt(sum((this_plate.ITs.temperature .^ 2) .* Pr) - avT^2);
                al = mean((this_plate.ITs.velocity .* this_plate.ITs.time)) ;
                al_mm = mean((this_plate.ITs.velocity_mm .* this_plate.ITs.time)) ;
                px2mm = (85 / (maxedge - minedge)); %85mm approx plate size / difference in edge clicks
                data_summary.plates.ITs_av_temp = [data_summary.plates.ITs_av_temp; avT];    
                data_summary.plates.ITs_av_temp_pos = [data_summary.plates.ITs_av_temp_pos; avT];    
                data_summary.plates.ITs_temp_band_width = [data_summary.plates.ITs_temp_band_width; BW]; 
                data_summary.plates.ITs_avg_length = [data_summary.plates.ITs_avg_length; al]; 
                data_summary.plates.ITs_avg_length_mm = [data_summary.plates.ITs_avg_length_mm; al_mm]; 
                data_summary.plates.px2mm = [data_summary.plates.px2mm; px2mm]; 
                data_summary.plates.ITs_total_time = [data_summary.plates.ITs_total_time; sum(this_plate.ITs.time)]; 
                data_summary.plates.ITs_total_dist = [data_summary.plates.ITs_total_dist; sum(this_plate.ITs.dist)]; 
                data_summary.plates.ITs_total_dist_mm = [data_summary.plates.ITs_total_dist_mm; sum(this_plate.ITs.dist_mm)];
                data_summary.plates.ITs_N = [data_summary.plates.ITs_N; length(this_plate.ITs.time)];   
                data_summary.plates.ITs_frac_time = [data_summary.plates.ITs_frac_time; sum(this_plate.ITs.time) ./  this_plate.plates.runs_total_time_all]; 
                data_summary.plates.ITs_frac_dist = [data_summary.plates.ITs_frac_dist; sum(this_plate.ITs.dist) ./  ...
                                         (this_plate.plates.runs_total_time_all .* (cur_plt_distance/cur_plt_time)) ]; %approx
                data_summary.plates.ITs_frac_N = [data_summary.plates.ITs_frac_N; length(this_plate.ITs.time)./ this_plate.plates.runs_N_all];
                
            else % no ITs on previous plate
                data_summary.plates.ITs_av_temp = [data_summary.plates.ITs_av_temp; -1];    
                data_summary.plates.ITs_temp_band_width = [data_summary.plates.ITs_temp_band_width; -1];
                data_summary.plates.ITs_avg_length = [data_summary.plates.ITs_avg_length; -1];
                data_summary.plates.ITs_avg_length_mm = [data_summary.plates.ITs_avg_length_mm; -1];
                data_summary.plates.px2mm = [data_summary.plates.px2mm; -1];
                data_summary.plates.ITs_total_time = [data_summary.plates.ITs_total_time; 0];   
                data_summary.plates.ITs_total_dist = [data_summary.plates.ITs_total_dist; 0];
                data_summary.plates.ITs_total_dist_mm = [data_summary.plates.ITs_total_dist_mm; 0];
                data_summary.plates.ITs_N = [data_summary.plates.ITs_N; 0];   
                data_summary.plates.ITs_frac_time = [data_summary.plates.ITs_frac_time; 0]; 
                data_summary.plates.ITs_frac_dist = [data_summary.plates.ITs_frac_dist; 0]; %approx
                data_summary.plates.ITs_frac_N = [data_summary.plates.ITs_frac_N; 0];
                
            end; % if length(this_plate.ITs.time) > 0
            
            this_plate.plates.turns_num_in_intrvl = [];
            for intrvl = 1:NUM_INTERVALS_FOR_ARS
                t1 = 300*(intrvl-1);
                t2 = 300*(intrvl);
                indx = find(this_plate.turns.t > t1 & this_plate.turns.t <= t2);
                this_plate.plates.turns_num_in_intrvl = [this_plate.plates.turns_num_in_intrvl; length(indx)];
            end;
            this_plate.plates.total_turns = sum(this_plate.plates.turns_num_in_intrvl);
            data_summary.plates.turns_num_in_intrvl = [data_summary.plates.turns_num_in_intrvl, this_plate.plates.turns_num_in_intrvl];
            data_summary.plates.total_turns = [data_summary.plates.total_turns; this_plate.plates.total_turns];
      	else % last plate was bad
            display(['Bad plate: ' current_plate]);
            % Runs    
            data_summary.runs.time = [data_summary.runs.time, {[]}]; 
            data_summary.runs.dx = [data_summary.runs.dx, {[]}]; 
        	data_summary.runs.angle = [data_summary.runs.angle, {[]}];
        	data_summary.runs.velocity = [data_summary.runs.velocity, {[]}];
            data_summary.runs.velocity_mm = [data_summary.runs.velocity_mm, {[]}];
            % ITs                    
          	data_summary.ITs.time = [data_summary.ITs.time, {[]}];       
          	data_summary.ITs.temperature = [data_summary.ITs.temperature, {[]}]; 
           	data_summary.ITs.velocity = [data_summary.ITs.velocity, {[]}];
            data_summary.ITs.velocity_mm = [data_summary.ITs.velocity_mm, {[]}];
            data_summary.ITs.X_pos = [data_summary.ITs.X_pos, {[]}];  
            data_summary.ITs.start_end = [data_summary.ITs.start_end, {[]}];  
            % Turns
          	data_summary.turns.t = [data_summary.turns.t, {[]}]; 
            data_summary.turns.x = [data_summary.turns.x, {[]}];   
            data_summary.turns.y = [data_summary.turns.y, {[]}];   
           	data_summary.turns.from_angle = [data_summary.turns.from_angle, {[]}]; 
          	data_summary.turns.to_angle = [data_summary.turns.to_angle, {[]}];
            % Plates (general, runs, ITs, turns)     
            data_summary.plates.plate_id = [data_summary.plates.plate_id, {current_plate}];
            data_summary.plates.av_velocity = [data_summary.plates.av_velocity; -999];
            data_summary.plates.av_velocity_mm = [data_summary.plates.av_velocity_mm; -999];
                                    %don't add bad plates to count
                                    
        	data_summary.plates.runs_total_time_dn = [data_summary.plates.runs_total_time_dn; -999];
        	data_summary.plates.runs_total_time_vr = [data_summary.plates.runs_total_time_vr; -999];
        	data_summary.plates.runs_total_time_up = [data_summary.plates.runs_total_time_up; -999];
            data_summary.plates.runs_total_time_all = [data_summary.plates.runs_total_time_all; -999];
        	data_summary.plates.runs_total_time_hr = [data_summary.plates.runs_total_time_hr; -999];
            
            data_summary.plates.runs_total_dx_dn = [data_summary.plates.runs_total_dx_dn; -999];
            data_summary.plates.runs_total_dx_vr = [data_summary.plates.runs_total_dx_vr; -999];
            data_summary.plates.runs_total_dx_up = [data_summary.plates.runs_total_dx_up; -999];

        	data_summary.plates.runs_N_dn = [data_summary.plates.runs_N_dn; -999];
        	data_summary.plates.runs_N_vr = [data_summary.plates.runs_N_vr; -999];
        	data_summary.plates.runs_N_up = [data_summary.plates.runs_N_up; -999];
            data_summary.plates.runs_N_all = [data_summary.plates.runs_N_all; -999];
        	data_summary.plates.runs_N_hr = [data_summary.plates.runs_N_hr; -999];

          	data_summary.plates.I_cr = [data_summary.plates.I_cr; -999]; 

            data_summary.plates.ITs_av_temp = [data_summary.plates.ITs_av_temp; -999];    
            data_summary.plates.ITs_temp_band_width = [data_summary.plates.ITs_temp_band_width; -999]; 
            data_summary.plates.ITs_avg_length = [  data_summary.plates.ITs_avg_length; -999];
            data_summary.plates.ITs_avg_length_mm = [  data_summary.plates.ITs_avg_length_mm; -999];
            data_summary.plates.px2mm = [  data_summary.plates.px2mm; -999];
            data_summary.plates.ITs_total_time = [data_summary.plates.ITs_total_time; -999]; 
            data_summary.plates.ITs_total_dist = [data_summary.plates.ITs_total_dist; -999];
            data_summary.plates.ITs_total_dist_mm = [data_summary.plates.ITs_total_dist_mm; -999];
            data_summary.plates.ITs_N = [data_summary.plates.ITs_N; -999];   
            
            data_summary.plates.turns_num_in_intrvl = [data_summary.plates.turns_num_in_intrvl, zeros(NUM_INTERVALS_FOR_ARS,1)-999];            
            data_summary.plates.total_turns = [data_summary.plates.total_turns, -999];            

        end; % last plate was bad
    end; % if rec_new_plate || rec_last_plate
    
    % organize data for current plate: 
    %---------------------------------
    if ~rec_last_plate
        maxedge = tr(k).max_edge;
        minedge = tr(k).min_edge;
        px2mm = (85 / (maxedge - minedge)); 
        % tr(k) is the current track - typically contains several runs and
        %                              turns, perhaps a few ITs. 
        current_plate = tr(k).file_name_id; % plate file-name prefix          
        x = tr(k).x; % x-axis position for each frame in the track
        y = tr(k).y; % y-axis position for each frame in the track
        t = tr(k).t; % absolute time (sec) for each frame in the track, e.g., t=[17,19,21,...] 
        r_time = tr(k).run_time; % duration in seconds of each run in the track
        r_dx = tr(k).run_dx;
        r_ang = tr(k).run_angle; % angle (rad) for each run in the track
        r_v = tr(k).run_velocities; % velocity (pxl/sec) during each run in the track
        r_v_mm = (tr(k).run_velocities)*px2mm; % velocity (mm/sec) during each run in the track
        IT_time = tr(k).IT_time; % duration (sec) for each IT in the track
        IT_dist = tr(k).IT_dist; % distance (px) for each IT in the track
        IT_dist_mm = (tr(k).IT_dist)*px2mm; % distance (mm) for each IT in the track
        IT_T = tr(k).IT_temperatures; % temperature (deg C) for each IT in the track 
        IT_X = tr(k).IT_X_postions;
        IT_v = tr(k).IT_velocities; % velocity (pxl/sec) during each IT in the track
        IT_v_mm = (tr(k).IT_velocities)*px2mm; % velocity (mm/sec) during each IT in the track
        maxT = tr(k).max_T;
        minT = tr(k).min_T;
        IT_ix =  tr(k).t(tr(k).IT_indx);
        if(numel(IT_ix) == 2)
            IT_ix = IT_ix';
        end
        maxedge = tr(k).max_edge;
        minedge = tr(k).min_edge;
        reversed = tr(k).reversed;
      	% Exploit symmetry to calculate if direction of run 
        % is up/down/perpendicular to the gradient
     	direction = r_ang; 
     	direction(direction>pi) = 2*pi-direction(direction>pi); 
        direction = abs(direction); 
      	% direction only has positive values between 0 and pi
            
        if rec_new_plate || k==1 % start organizing data for new plate
            if isempty(t) % initialize distance & time for new plate
                cur_plt_distance = 0;
                cur_plt_distance_mm = 0;
                cur_plt_time = 0;
            else
                cur_plt_distance = sum(sqrt(diff(x).^2 + diff(y).^2));
                cur_plt_distance_mm = (sum(sqrt(diff(x).^2 + diff(y).^2)))*px2mm;
                cur_plt_time = t(end)-t(1)+1;
            end;        
            % Runs
            this_plate.runs.time = r_time; 
            this_plate.runs.dx = r_dx; 
            this_plate.runs.angle = r_ang;
            this_plate.runs.velocity = r_v;
            this_plate.runs.velocity_mm = r_v_mm;
            % ITs
            this_plate.ITs.time = IT_time;  
            this_plate.ITs.dist = IT_dist;
            this_plate.ITs.dist_mm = IT_dist_mm;
            this_plate.ITs.temperature = IT_T; 
            this_plate.ITs.X_pos = IT_X;
            this_plate.ITs.start_end = IT_ix;
            this_plate.ITs.velocity = IT_v;                
            this_plate.ITs.velocity_mm = IT_v_mm;
            % Turns 
            this_plate.turns.t = t(tr(k).turn_indx);  
            this_plate.turns.x = x(tr(k).turn_indx);  
            this_plate.turns.y = y(tr(k).turn_indx);  

            if ~isempty(this_plate.turns.t)
                this_plate.turns.from_angle = r_ang(1:(end-1)); 
                this_plate.turns.to_angle = r_ang(2:end); 
            else
                this_plate.turns.from_angle = [];
                this_plate.turns.to_angle = [];
            end;           
            % Plates (runs, ITs)
            
            indx_dn = find(direction>(pi-ANG_FOR_RW) & direction<=pi);
            if isempty(indx_dn)
                this_plate.plates.runs_total_time_dn = 0;
                this_plate.plates.runs_total_dx_dn = 0;
                this_plate.plates.runs_N_dn = 0;
            else
                this_plate.plates.runs_total_time_dn = sum(r_time(indx_dn));
                this_plate.plates.runs_total_dx_dn = sum(r_dx(indx_dn));
                this_plate.plates.runs_N_dn = length(indx_dn);
            end;
            indx_vr = find(direction<(pi/2+ANG_FOR_RW) & direction>pi/2-ANG_FOR_RW);
            if isempty(indx_vr)
                this_plate.plates.runs_total_time_vr = 0;
                this_plate.plates.runs_total_dx_vr = 0;
                this_plate.plates.runs_N_vr = 0;
            else
                this_plate.plates.runs_total_time_vr = sum(r_time(indx_vr));
                this_plate.plates.runs_total_dx_vr = sum(r_dx(indx_vr));
                this_plate.plates.runs_N_vr = length(indx_vr);
            end;
            indx_up = find(direction<ANG_FOR_RW & direction>=0);
            if isempty(indx_up)
                this_plate.plates.runs_total_time_up = 0;
                this_plate.plates.runs_total_dx_up = 0;
                this_plate.plates.runs_N_up = 0;
            else
                this_plate.plates.runs_total_time_up = sum(r_time(indx_up));
                this_plate.plates.runs_total_dx_up = sum(r_dx(indx_up));
                this_plate.plates.runs_N_up = length(indx_up);
            end;

                this_plate.plates.runs_total_time_all = this_plate.plates.runs_total_time_up  +  this_plate.plates.runs_total_time_dn +  this_plate.plates.runs_total_time_vr;
                this_plate.plates.runs_total_time_hr = this_plate.plates.runs_total_time_up  +  this_plate.plates.runs_total_time_dn;
                this_plate.plates.runs_N_all = this_plate.plates.runs_N_up + this_plate.plates.runs_N_dn + this_plate.plates.runs_N_vr;
                this_plate.plates.runs_N_hr = this_plate.plates.runs_N_up + this_plate.plates.runs_N_dn;


            
            if isempty(IT_time)
                this_plate.plates.ITs_total_time = 0; 
                this_plate.plates.ITs_N = 0;         
                this_plate.plates.ITs_total_dist = 0;
                this_plate.plates.ITs_total_dist_mm = 0;
                
            else
                this_plate.plates.ITs_total_time = sum(IT_time); 
                this_plate.plates.ITs_N = length(IT_time);       
                this_plate.plates.ITs_total_dist = sum(IT_dist); 
                this_plate.plates.ITs_total_dist_mm = sum(IT_dist_mm);
            end;
                    
        else % continuing with current plate
                       
            if ~isempty(t) 
                cur_plt_distance = cur_plt_distance + sum(sqrt(diff(x).^2 + diff(y).^2));
                cur_plt_distance_mm = cur_plt_distance_mm + (sum(sqrt(diff(x).^2 + diff(y).^2)))*px2mm;
                cur_plt_time = cur_plt_time + t(end)-t(1)+1;
            end;        
            % Runs
            this_plate.runs.time = [this_plate.runs.time; r_time];   
            this_plate.runs.dx = [this_plate.runs.dx; r_dx];
            this_plate.runs.angle = [this_plate.runs.angle; r_ang];
            this_plate.runs.velocity = [this_plate.runs.velocity; r_v];
            this_plate.runs.velocity_mm = [this_plate.runs.velocity_mm; r_v_mm];
            % ITs
            this_plate.ITs.time = [this_plate.ITs.time; IT_time];      
            this_plate.ITs.dist = [this_plate.ITs.dist; IT_dist];
            this_plate.ITs.dist_mm = [this_plate.ITs.dist_mm; IT_dist_mm];
            this_plate.ITs.temperature = [this_plate.ITs.temperature; IT_T]; 
            this_plate.ITs.X_pos = [this_plate.ITs.X_pos; IT_X]; 
            this_plate.ITs.start_end = [this_plate.ITs.start_end; IT_ix]; 
            this_plate.ITs.velocity = [this_plate.ITs.velocity; IT_v]; 
            this_plate.ITs.velocity_mm = [this_plate.ITs.velocity_mm; IT_v_mm]; 
            % Turns 
            this_plate.turns.t = [this_plate.turns.t; t(tr(k).turn_indx)];   
            this_plate.turns.x = [this_plate.turns.x; x(tr(k).turn_indx)];   
            this_plate.turns.y = [this_plate.turns.y; y(tr(k).turn_indx)];   
            if ~isempty(this_plate.turns.t)
                this_plate.turns.from_angle = [this_plate.turns.from_angle; r_ang(1:(end-1))]; 
                this_plate.turns.to_angle = [this_plate.turns.to_angle; r_ang(2:end)]; 
            end;           
            % Plates (runs and ITs)
            indx_dn = find(direction>(pi-ANG_FOR_RW) & direction<=pi);
            if ~isempty(indx_dn)
                this_plate.plates.runs_total_time_dn = this_plate.plates.runs_total_time_dn+sum(r_time(indx_dn));
                this_plate.plates.runs_total_dx_dn = this_plate.plates.runs_total_dx_dn+sum(r_dx(indx_dn));
                this_plate.plates.runs_N_dn = this_plate.plates.runs_N_dn+length(indx_dn);
            end;
            indx_vr = find(direction<(pi/2+ANG_FOR_RW) & direction>pi/2-ANG_FOR_RW);
            if ~isempty(indx_vr)
                this_plate.plates.runs_total_time_vr = this_plate.plates.runs_total_time_vr+sum(r_time(indx_vr));
                this_plate.plates.runs_total_dx_vr = this_plate.plates.runs_total_dx_vr+sum(r_dx(indx_vr));
                this_plate.plates.runs_N_vr = this_plate.plates.runs_N_vr+length(indx_vr);
            end;
            indx_up = find(direction<ANG_FOR_RW & direction>=0);
            if ~isempty(indx_up)
                this_plate.plates.runs_total_time_up = this_plate.plates.runs_total_time_up+sum(r_time(indx_up));
                this_plate.plates.runs_total_dx_up = this_plate.plates.runs_total_dx_up+sum(r_dx(indx_up));
                this_plate.plates.runs_N_up = this_plate.plates.runs_N_up+length(indx_up);
            end;

            
                this_plate.plates.runs_total_time_all = this_plate.plates.runs_total_time_up  +  this_plate.plates.runs_total_time_dn +  this_plate.plates.runs_total_time_vr;
                this_plate.plates.runs_total_time_hr = this_plate.plates.runs_total_time_up  +  this_plate.plates.runs_total_time_dn;
                this_plate.plates.runs_N_all = this_plate.plates.runs_N_up + this_plate.plates.runs_N_dn + this_plate.plates.runs_N_vr;
                this_plate.plates.runs_N_hr = this_plate.plates.runs_N_up + this_plate.plates.runs_N_dn;


            
            

            if ~isempty(IT_time)
                this_plate.plates.ITs_total_time = this_plate.plates.ITs_total_time+sum(IT_time); 
                this_plate.plates.ITs_N = this_plate.plates.ITs_N+length(IT_time); 
                this_plate.plates.ITs_total_dist = this_plate.plates.ITs_total_dist+sum(IT_dist); 
                this_plate.plates.ITs_total_dist_mm = this_plate.plates.ITs_total_dist_mm+sum(IT_dist_mm); 
            end;
        end; % continuing with current plate
        
    end; % if ~rec_last_plate
    
end; % for k = 1:(length(tr)+1)ma

% Average data from all plates
%-----------------------------
N = data_summary.plates.num_of_plates;
sqrtN = sqrt(N);

m = mean(data_summary.plates.av_velocity);
e = std(data_summary.plates.av_velocity)/sqrtN; 
data_summary.ALL.all_mean_velocity = [m, e]; 

m = mean(data_summary.plates.av_velocity_mm);
e = std(data_summary.plates.av_velocity_mm)/sqrtN; 
data_summary.ALL.all_mean_velocity_mm = [m, e]; 

plts = [];
plts_mm = [];

for k=1:N
    tmp = data_summary.runs.velocity{k};
    tmp_mm = data_summary.runs.velocity_mm{k};
    if ~isempty(tmp); % not bad plate
        plts = [plts mean(tmp)]; % mean of single plate
        plts_mm = [plts_mm mean(tmp_mm)]; % mean of single plate
    end;
end;
m = mean(plts);
e = std(plts)/sqrtN; 
data_summary.ALL.run_mean_velocity = [m, e]; 
m = mean(plts_mm);
e = std(plts_mm)/sqrtN; 
data_summary.ALL.run_mean_velocity_mm = [m, e]; 

m = mean(data_summary.plates.runs_total_time_dn);
e = std(data_summary.plates.runs_total_time_dn)/sqrtN;
data_summary.ALL.run_dn_total_duration = [m, e];
m = mean(data_summary.plates.runs_total_time_vr);
e = std(data_summary.plates.runs_total_time_vr)/sqrtN;
data_summary.ALL.run_vr_total_duration = [m, e]; 
m = mean(data_summary.plates.runs_total_time_up);
e = std(data_summary.plates.runs_total_time_up)/sqrtN;
data_summary.ALL.run_up_total_duration = [m, e]; 
m = mean(data_summary.plates.runs_total_time_all);
e = std(data_summary.plates.runs_total_time_all)/sqrtN;
data_summary.ALL.run_all_total_duration = [m, e];
m = mean(data_summary.plates.runs_total_time_hr);
e = std(data_summary.plates.runs_total_time_hr)/sqrtN;
data_summary.ALL.run_hr_total_duration = [m, e];

m = mean(data_summary.plates.runs_N_dn);
e = std(data_summary.plates.runs_N_dn)/sqrtN;
data_summary.ALL.run_dn_num = [m, e];            
m = mean(data_summary.plates.runs_N_vr);
e = std(data_summary.plates.runs_N_vr)/sqrtN;
data_summary.ALL.run_vr_num = [m, e];            
m = mean(data_summary.plates.runs_N_up);
e = std(data_summary.plates.runs_N_up)/sqrtN;
data_summary.ALL.run_up_num = [m, e];

m = mean(data_summary.plates.runs_N_all);
e = std(data_summary.plates.runs_N_all)/sqrtN;
data_summary.ALL.run_all_num = [m, e];
m = mean(data_summary.plates.runs_N_hr);
e = std(data_summary.plates.runs_N_hr)/sqrtN;
data_summary.ALL.run_hr_num = [m, e];


m = mean(data_summary.plates.runs_total_time_dn./data_summary.plates.runs_N_dn);
e = std(data_summary.plates.runs_total_time_dn./data_summary.plates.runs_N_dn)/sqrtN;
data_summary.ALL.run_dn_mean_duration = [m, e]; 
m = mean(data_summary.plates.runs_total_time_vr./data_summary.plates.runs_N_vr);
e = std(data_summary.plates.runs_total_time_vr./data_summary.plates.runs_N_vr)/sqrtN;
data_summary.ALL.run_vr_mean_duration = [m, e]; 
m = mean(data_summary.plates.runs_total_time_up./data_summary.plates.runs_N_up);
e = std(data_summary.plates.runs_total_time_up./data_summary.plates.runs_N_up)/sqrtN;
data_summary.ALL.run_up_mean_duration = [m, e]; 
m = mean(data_summary.plates.runs_total_time_all./data_summary.plates.runs_N_all);
e = std(data_summary.plates.runs_total_time_all./data_summary.plates.runs_N_all)/sqrtN;
data_summary.ALL.run_all_mean_duration = [m, e]; 
m = mean(data_summary.plates.runs_total_time_hr./data_summary.plates.runs_N_hr);
e = std(data_summary.plates.runs_total_time_hr./data_summary.plates.runs_N_hr)/sqrtN;
data_summary.ALL.run_hr_mean_duration = [m, e]; 

m = mean(data_summary.plates.I_cr);
e = std(data_summary.plates.I_cr)/sqrtN;
data_summary.ALL.I_cr = [m, e];  

m = mean(data_summary.plates.I_cr_dx);
e = std(data_summary.plates.I_cr_dx)/sqrtN;
data_summary.ALL.I_cr_dx = [m, e]; 
 

plts = [];
for k=1:N
    tmp = data_summary.ITs.velocity{k};
    temp_mm = data_summary.ITs.velocity_mm{k};
    if ~isempty(tmp); % not bad plate
        plts = [plts mean(tmp)]; % mean of single plate
        plts_mm = [plts_mm mean(tmp_mm)]; % mean of single plate
    end;
end;
m = mean(plts);
e = std(plts)/sqrt(length(plts)); % perhaps not every plate had ITs
data_summary.ALL.IT_mean_velocity = [m, e];
m = mean(plts_mm);
e = std(plts_mm)/sqrt(length(plts_mm)); % perhaps not every plate had ITs
data_summary.ALL.IT_mean_velocity_mm = [m, e];

tmp1 = data_summary.plates.ITs_total_time;
tmp1 = tmp1(tmp1>0); % perhaps not every plate had ITs
m = mean(tmp1); 
e = std(tmp1)/sqrt(length(tmp1));
data_summary.ALL.IT_total_duration = [m, e];  

tmp2 = data_summary.plates.ITs_N;
tmp2 = tmp2(tmp2>0); % perhaps not every plate had ITs
m = mean(tmp2); 
e = std(tmp2)/sqrt(length(tmp2));
data_summary.ALL.IT_num = [m, e];   


tmp2 = data_summary.plates.ITs_frac_time;
tmp2 = tmp2(tmp2>0); % perhaps not every plate had ITs
m = mean(tmp2); 
e = std(tmp2)/sqrt(length(tmp2));
data_summary.ALL.ITs_frac_time = [m, e];

tmp2 = data_summary.plates.ITs_frac_dist;
tmp2 = tmp2(tmp2>0); % perhaps not every plate had ITs
m = mean(tmp2); 
e = std(tmp2)/sqrt(length(tmp2));
data_summary.ALL.ITs_frac_dist = [m, e];


tmp2 = data_summary.plates.ITs_frac_N;
tmp2 = tmp2(tmp2>0); % perhaps not every plate had ITs
m = mean(tmp2); 
e = std(tmp2)/sqrt(length(tmp2));
data_summary.ALL.ITs_frac_N = [m, e];

tmp1 = data_summary.plates.ITs_total_time;
tmp1 = tmp1(tmp1>0); % perhaps not every plate had ITs
tmp2 = data_summary.plates.ITs_N;
tmp2 = tmp2(tmp2>0); % perhaps not every plate had ITs
%fix
m = mean(tmp1./tmp2); 
e = std(tmp1./tmp2)/sqrt(length(tmp1));
data_summary.ALL.IT_mean_duration = [m, e]; 

tmp = data_summary.plates.ITs_av_temp;
tmp = tmp(tmp>0); % perhaps not every plate had ITs
m = mean(tmp); 
e = std(tmp)/sqrt(length(tmp));
data_summary.ALL.IT_mean_temperature = [m, e];

tmp = data_summary.plates.ITs_temp_band_width;
tmp = tmp(tmp>0); % perhaps not every plate had ITs
m = mean(tmp); 
e = std(tmp)/sqrt(length(tmp));
data_summary.ALL.IT_band_width = [m, e]; 

tmp = data_summary.plates.ITs_avg_length;
tmp = tmp(tmp>0); % perhaps not every plate had ITs
m = mean(tmp); 
e = std(tmp)/sqrt(length(tmp));
data_summary.ALL.ITs_avg_length = [m, e]; 

tmp_mm = data_summary.plates.ITs_avg_length_mm;
tmp_mm = tmp_mm(tmp_mm>0); % perhaps not every plate had ITs
m = mean(tmp_mm); 
e = std(tmp_mm)/sqrt(length(tmp_mm));
data_summary.ALL.ITs_avg_length_mm = [m, e]; 


tmp = data_summary.plates.px2mm;
m = mean(data_summary.plates.px2mm);
e = std(tmp)/sqrt(length(tmp));
data_summary.ALL.px2mm = [m, e]; 

m = mean(data_summary.plates.turns_num_in_intrvl,2); % average each row (5 minute time-interval) over all plates
e = std(data_summary.plates.turns_num_in_intrvl,0,2)/sqrtN; % standard error for each row (5 minute time-interval) 
data_summary.ALL.turns_num_in_intrvl = [m, e]; % each row corresponds to a 5 min interval

m = mean(data_summary.plates.total_turns); % average each row (5 minute time-interval) over all plates
e = std(data_summary.plates.total_turns)/sqrtN; 
data_summary.ALL.total_turns = [m, e]; 
%%%%%%%%%%%%%%%%%%%%%%
%%% Make Log file: %%%
%%%%%%%%%%%%%%%%%%%%%%

% text file
dir_name = data_summary.path; 
tmp = strfind(dir_name,filesep);
sv_nm = ['LG_' dir_name((tmp(end)+1):end) '_'];
tmp = datestr(now,31);
tmp = ['D' tmp(1:end-3)]; % get rid of seconds
tmp = regexprep(tmp,{':'},'_');
tmp = regexprep(tmp,{' '},'_H');
sv_nm = [sv_nm tmp '.txt'];
data_sv_nm = regexprep(sv_nm,'LG_','DS_');
data_sv_nm = regexprep(data_sv_nm,'.txt','.mat');

fh = fopen([dir_name filesep sv_nm],'w');
fprintf(fh,' \r\n Log file name  : %s\r\n',sv_nm);
fprintf(fh,' Data file name : %s\r\n',data_sv_nm);
fprintf(fh,' Path: %s\r\n',data_summary.path);
fprintf(fh,' Mean Relative Error: %s\r\n',num2str(data_summary.mean_relative_error));
fprintf(fh,' \r\n');
fprintf(fh,'Data per plate: \r\n--------------- \r\n');
fns = fieldnames(data_summary.plates);
for k = 1:length(fns)
    if(isequal(fns{k},'max_T')||isequal(fns{k},'min_T')||isequal(fns{k},'max_edge')||isequal(fns{k},'reversed')||isequal(fns{k},'min_edge')||isequal(fns{k},'units')||isequal(fns{k},'ITs_av_temp_pos'))
       continue; 
    end
    nxt_field = data_summary.plates.(fns{k}); 
    if isequal(fns{k},'turns_num_in_intrvl')
        for j = 1:size(nxt_field,1)
            m1 = (j-1)*5+1;
            m2 = j*5; 
            [s, e] = sprintf('%s min %s-%s :',fns{k},num2str(m1),num2str(m2));
            fprintf(fh,'%-32s %s \r\n',s,num2str(nxt_field(j,:)));
        end; 
        
    else
        if isnumeric(nxt_field)
            data_str = num2str(nxt_field'); % make it a row
        elseif iscell(nxt_field) 
            data_str = [];
            for j = 1:length(nxt_field)
                data_str = [data_str, ' ', nxt_field{j} ' '];
            end;
        end;
        fprintf(fh,'%-19s: %s',fns{k},data_str);
        fprintf(fh,' %s \r\n',getfield(data_summary.plates.units,fns{k}));
    end;
end;
fprintf(fh,' \r\n');
fprintf(fh,'Averaged over all plates [mean, standard-error]: \r\n------------------------------------------------ \r\n');
fns = fieldnames(data_summary.ALL);
for k = 1:length(fns)
    if(isequal(fns{k},'max_T')||isequal(fns{k},'min_T')||isequal(fns{k},'max_edge')||isequal(fns{k},'reversed')||isequal(fns{k},'min_edge')||isequal(fns{k},'units')||isequal(fns{k},'ITs_av_temp_pos'))
       continue; 
    end
    nxt_field = data_summary.ALL.(fns{k}); 
    
    if isequal(fns{k},'turns_num_in_intrvl')
    	for j = 1:size(nxt_field,1)
            m1 = (j-1)*5+1;
            m2 = j*5; 
            [s, e] = sprintf('%s min %s-%s :',fns{k},num2str(m1),num2str(m2));
            fprintf(fh,'%-32s %s \r\n',s,num2str(nxt_field(j,:)));
        end;
    else
        data_str = num2str(nxt_field); % make it a row
        fprintf(fh,'%-21s : %s',fns{k},data_str);
        fprintf(fh,' %s \r\n',getfield(data_summary.ALL.units,fns{k}));
    end;
end;
fprintf(fh,' \r\n');
fclose(fh); 

%s = textread([dir_name '\' sv_nm], '%s', 'whitespace',''); % preserve white-spaces
%display(s{1}); 
data_summary.tr = tr;
% data file 
fh = fopen([dir_name filesep data_sv_nm],'w');
save([dir_name filesep data_sv_nm], 'data_summary','-mat');
fclose(fh); 

type([dir_name filesep sv_nm])

return; 



