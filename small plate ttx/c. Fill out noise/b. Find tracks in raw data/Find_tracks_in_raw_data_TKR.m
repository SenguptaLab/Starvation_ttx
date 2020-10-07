close all;
clear all;
clc;

%%% for IT_index calculate by hand:
%%% Baseline = (N2 @ 0.7 deg/cm) - (N2 @ 0 deg/cm)
%%% Index = [(mutant @ 0.7 deg/cm) - (mutant @ 0 deg/cm)] / baseline
%%% (calculate for average total time spents tracking per plate, 
%%%                average number of tracks per plate, 
%%%                avrage dura tion of ITs per plate (= ratio of total time / number of tracks)

%%% IMPORTANT - OUR DATA IS SOMETIMES SAMPLED AT 1 FPS %%%
%%%            DAMON'S IS USUALLY SAMPLED AT 0.5 FPS   %%%
global FPS
global SKIP_EVEN_LINES
newcam=2;
%FPS = 0.5; % rate of image capturing in Frames Per Second 
FPS = 1.0; % 1 for Sara/Matt data; 0.5 for Alex-Damon data
SKIP_EVEN_LINES = 0; % set to 1 for emulating 1/2 FPS on 1 FPS centers file

%%%%%%%% ONLY ANALYZE T_START < time < T_END %%%%%%%%
%%%  300 sec =  5 minutes                         %%%
%%%  600 sec = 10 minutes                         %%%
%%%  900 sec = 15 minutes                         %%%
%%% 1200 sec = 20 minutes                         %%%
%%% 1500 sec = 25 minutes                         %%%
%%% 1800 sec = 30 minutes                         %%%
%%% 2100 sec = 35 minutes                         %%%
%%% 2400 sec = 40 minutes                         %%%
%%% 2700 sec = 45 minutes                         %%%
%%% 3000 sec = 50 minutes                         %%%
FRM_START = 0*FPS; % number is in sec           %%%
%FRM_START = 1*FPS; % number is in sec           %%%

                     % number * FPS = # of frames %%%
%FRM_END = 2099*FPS;  % Only analyze cryophilic    %%%
FRM_END = 3813*FPS;  % Only analyze cryophilic    %%%
%FRM_END = 900*FPS;  % Only analyze cryophilic    %%%

                     % migration time             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_mode = 'DB_global'; % Use David's flagturnsGlobal4 
%run_mode = 'DB_local'; % Use David's flagturnsLlocal4 

%%% Checking stuff with Damon %%%
%%%% THESE SHOULD ONLY BE CHANGED FOR DEBUGGING %%%
QUICK_RUN_2_N2_PLATES = 0; % Should be 0 when real data is analyzed
%run_mode = 'DBGpi/4'; % Use David's flagturnsGlobal4 and set ALPHA=pi/4
%run_mode = 'DBLpi/4'; % Use David's flagturnsLocal4 and set ALPHA=pi/4
%run_mode = 'AlDa'; % Use Alex's and Damon's flagturns
%run_mode = 'DB_AlDaCutOffs'; % Use David's Local turnd code with cutoffs inspired by Alex-Damon code
%run_mode = 'DB_Global_CutOffs'; % Use David's Global turnd code with cutoffs inspired by Alex-Damon code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%basedir = 'C:\Documents and Settings\smwasser\My Documents\Work\';
%basedir = 'Y:\David Test\';
%restoredefaultpath; % restore default path
basedir =  fileparts(pwd);

datadir = [basedir 'z. data/'];
pth = genpath(basedir); % returns a path string formed by recursively 
                        % adding all the directories below basedir
addpath(pth);
t0 = clock; 

%%%%%%%%% Parameters for all scripts  (GLOBAL) %%%%%%%%%%%%%%

% readallcenters.m: (parameters for the track.m function)
%-------------------
global MAX_DIST 
global PARAM_MEM
global PARAM_DIM
global PARAM_GOOD
global PARAM_QUIET

MAX_DIST = 5/FPS; % max no. of pixels worms can move between consecutive frames
PARAM_MEM = 1; % max size of "lost frames gap" allowed in a single particle track
PARAM_DIM = 2; % dimensions of tracking problem (plate --> 2D)
PARAM_GOOD = 4*FPS; % minimal no. of frames in "legal" track - shorter are discarded

PARAM_QUIET = 1; % allow track.m to display text comments (set to 0 for comments, 1 for quiet)
%--------------------

% putinfields.m:
%---------------
global TOOSHORT
global MAXGAP
global DEFAULT_MOVEBY

TOOSHORT = PARAM_GOOD; % minimal number of frames in "legal" track
MAXGAP = 1; % max frame gap allowed in single particle track 
            % (more than that gets split to two tracks)
            % THIS SHOULD NEVER BE LARGER THAN 1
DEFAULT_MOVEBY = 5*newcam; % min movement in pixels by worm in either x- or 
                    % y-direction in single track 
                    % (less is considered noise / dead worm)
%---------------

%flagturns4.m: 
%-------------
global ALPHA   % minimal changle in angle that counts for a turn

global SEG_D   % distance in pixels that counts as a discrete "step"
global SEG_T1  % minimal time in seconds for a segment ro be counted
global SEG_T2  % maximal time in seconds for a segment ro be counted
global N_SEG   % minimal number of segments for "run"  
global MAXDANG_IT % max deviation from pi/2 (or -pi/2) which is still considered IT (radians)
global MINLENGTH_IT % min time of track to be considered IT (in seconds)

ALPHA  = pi/6.5;  % minimal changle in angle that counts for a turn

%SEG_D  =  4;    % distance in pixels that counts as a discrete "segment"
SEG_D  =  4*newcam;    % distance in pixels that counts as a discrete "segment"
SEG_T1 =  3;    % minimal time in seconds for a segment to be counted
SEG_T2 = 10;    % maximal time in seconds for a segment to be counted
N_SEG  =  2;    % minimal number of segments for "run"  
MAXDANG_IT = pi/30; % max deviation from pi/2 (or -pi/2) which is still considered IT (radians)
%%% MThreshold time in seconds used to determine if a %%%
%%% verticle track is an IT. Default value = 20 sec   %%%
%%% (20sec ~ 4 body lengths ~ fairly short tracks     %%%
MINLENGTH_IT = 3.6;    % min length of track to be considered IT (in mm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------

% make_data_summary4.m and make_plot functions
%---------------------------------------------
global TIME_BIN
global ANG_BIN
global MAX_RUN_DURATION_FOR_RW
global MIN_RUN_DURATION_FOR_RW
global ANG_FOR_RW

TIME_BIN = 8; % bin width for "runs" times in runs histogams (sec)
ANG_BIN = pi/8; % bin width for angles in run time vs. angle plot (rad)
MAX_RUN_DURATION_FOR_RW = 600; % seconds
MIN_RUN_DURATION_FOR_RW = 2; % seconds
ANG_FOR_RW = pi/4;
%-------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Convert centers file to particle tracks data (and save data in file) %%% 
if QUICK_RUN_2_N2_PLATES
    dir_name = 'C:\Documents and Settings\wormly\Desktop\David Test\z. data\Sara_Old_N2_Cryo_2_plates';
else
    dir_name = uigetdir(datadir); 
end;

tmp1 = readallcenters4(dir_name,FRM_START,FRM_END); % reads center files, tracks particles, 
                              % arranges data in tr structure format    
                              
                              

                              
if isequal(run_mode,'DB_global')
    tr = flagturnsGlobal4(tmp1); % flag turns, run durations, run lengths etc.
    clear('tmp1');
%    save_tr4(dir_name,'OdbG_',tr);       % tr has all fields
elseif isequal(run_mode,'DB_local')
    tr = flagturnsLocal4(tmp1); % flag turns, run durations, run lengths etc.
    clear('tmp1');
%    save_tr4(dir_name,'OdbL_',tr);       % tr has all fields
elseif isequal(run_mode,'DBGpi/4')
    ALPHA = pi/4; 
    tr = flagturnsGlobal4(tmp1); % angle for turn is pi/4
    clear('tmp1');
%    save_tr4(dir_name,'OdbGp4_',tr);       % tr has all fields
elseif isequal(run_mode,'DBLpi/4')
    ALPHA = pi/4; 
    tr = flagturnsLocal4(tmp1); % angle for turn is pi/4
    clear('tmp1');
%    save_tr4(dir_name,'OdbLp4_',tr);       % tr has all fields
elseif isequal(run_mode,'AlDa')
    tr = AlexDamon_algorithm_flagturns(tmp1);
    clear('tmp1');
%    save_tr4(dir_name,'Oalda_',tr);       % tr has all fields
elseif isequal(run_mode,'DB_AlDaCutOffs')
    ALPHA = pi/4; 
    tr = flagturnsLocal_AlDaCutoffs4(tmp1); % angle for turn is pi/4
    clear('tmp1');
%    save_tr4(dir_name,'AlDaCutoffs_',tr);       % tr has all fields
elseif isequal(run_mode,'DB_Global_CutOffs')
    ALPHA = pi/6; 
    tr = flagturnsGlobal_Cutoffs4(tmp1); % angle for turn is pi/4
    clear('tmp1');
%    save_tr4(dir_name,'AlDaCutoffs_',tr);       % tr has all fields
else
    display(['Error: invalid run_mode - ' run_mode]);
    return;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% plot example tracks and display run time %%% 
%Plot_Tracks_and_turns(tr);
display(['Run time: ' num2str(round(etime(clock,t0)/6)/10) ' minutes.']);  % time in minutes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return;

