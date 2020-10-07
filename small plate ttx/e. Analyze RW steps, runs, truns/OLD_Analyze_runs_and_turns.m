% Run this script to Analyze runs and turns.
% It sets the valuse of global parameters necessary for the analysis. 

%Global parameters for runlengthanal4.m: 
global TIME_BIN
global ANG_BIN_1
global MAX_RUN_DURATION_FOR_RW
global MIN_RUN_DURATION_FOR_RW
global ANG_FOR_RW

TIME_BIN = 8; % bin-width for "runs" duration times in seconds
ANG_BIN_1 = pi/8;
MAX_RUN_DURATION_FOR_RW = 600; % seconds
MIN_RUN_DURATION_FOR_RW = 2; % seconds
ANG_FOR_RW = pi/4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PLATE_LIST = []; % list of plate numbers to analyze ([] means analyze all)
t0 = clock; 

basedir = 'Y:\David Test\z. data\';
% Only load tr from file if not in memory! 
if isempty(whos('tr')) % This means that tr is not in the memory and needs uploading
    display(['Error: ''tr'' structure not found. ']);
    display(['Please run ''Find_tracks_in_raw_data.m'' and then retry.']);
    return;
    %%% Older code designed to get tr from a mat file where it was %%%
    %%% saved. Newer code does not save tr (too big) but only a    %%% 
    %%% data-summary structure.                                    %%%  
    %[fn, pth] = uigetfile([basedir '*.mat'],'Select tracks data file:');
    %tmp1 = load([pth '\' fn]);       % should give a structure: tmp1.tr 
    %tmp2 = fieldnames(tmp1);         % should be only one field, called "tr"
    %tr = tmp1.(cell2mat(tmp2(1)));   % tmp1.tr contains the tracks data
    %clear('tmp1','tmp2');            % save memory, tmo1 and tr might be very large
else % tr is in memory already, just choose path
    pth = tr(1).path;
end;

%%% Determine which plates %%%
plates_indx = [];
plate_nums = [];
for k=1:length(tr)
    cur_id = tr(k).plate_id; 
	flag = ismember(cur_id,PLATE_LIST) || isempty(PLATE_LIST);
    if  flag % find index in tr2 of tracks from listed plates (or all if list is empty)
    	plates_indx = [plates_indx, k];
        plate_nums = [plate_nums cur_id];
	end;
end;
plate_nums = unique(plate_nums);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% plot times of tracks %%%%%%%%%%%%%
%for i = 1:length(plate_nums)
%    figure(2000+i);
%    clf;
%end; 
%for i_indx = 1:length(plates_indx)
%    i = plates_indx(i_indx);
%    t = tr(i).t;
%    figure(2000+tr(i).plate_id);
%    hold on;
%    plot([t(1), t(end)],[i,i]);
%    title(['Cleaned (real) tracks found for plate no. ' num2str(tr(i).plate_id)]);
%    hold off;
%end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

Icr_by_plate = [];
Ns_by_plate = [];
Durs_by_plate = [];
Vs_by_plate = [];
IT_t = [];
IT_N = [];
IT_sum_t = [];
IT_time_fracs = [];
IT_dist_fracs = [];
IT_center_T = [];
IT_width_T = [];
for p = 1:length(plate_nums)
    [Icr,Ns,Durs, Vs] = runlengthanal4(tr,[p]);
    Icr_by_plate = [Icr_by_plate; Icr];
    Ns_by_plate = [Ns_by_plate; Ns]; % each Ns is row of 3
    Durs_by_plate = [Durs_by_plate; Durs]; % each Durs is row of 3
    Vs_by_plate = [Vs_by_plate; Vs]; % each Vs is row of 3
    
    % IT data for current plate
    all_plt_indx = [];
    for k = 1:length(tr)
        all_plt_indx = [all_plt_indx; tr(k).plate_id];
    end;
    indx = find(all_plt_indx == plate_nums(p));
    cur_tr = tr(indx);
    all_IT_N = 0;
    all_IT_t = [];
    all_IT_sum_t = 0;
    all_IT_dist = [];
    all_IT_Temp = [];
    all_run_t = [];
    all_run_dist = [];
    for k = 1:length(cur_tr) % run over all tracks in current plate
        all_IT_t = [all_IT_t; cur_tr(k).IT_time];
        all_IT_N = all_IT_N+length(cur_tr(k).IT_time);
        all_IT_sum_t = all_IT_sum_t+sum(cur_tr(k).IT_time);
        all_IT_dist = [all_IT_dist; cur_tr(k).IT_dist];
        all_IT_Temp = [all_IT_Temp; cur_tr(k).IT_temperatures];
        all_run_t = [all_run_t; cur_tr(k).run_time];
        all_run_dist = [all_run_dist; cur_tr(k).run_dist];
    end;  
    all_run_t = all_run_t(all_run_t > MINTIME_IT); % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cur_file_num = cur_tr(1).file_name_id; 
    all_IT_time_frac = sum(all_IT_t) / sum(all_run_t); % fraction of IT time
    all_IT_dist_frac = sum(all_IT_dist) / sum(all_run_dist); % fraction of IT distance
    all_IT_av_T = sum(all_IT_Temp .* all_IT_t) / sum(all_IT_t); % weighted (by duration) average Temparature of all ITs
    sncd_moment = sum((all_IT_Temp.^2) .* all_IT_t) / sum(all_IT_t); % second moment of IT Temperatures, weighted by duration 
    all_IT_std_T = sqrt(sncd_moment - all_IT_av_T^2); % width IT band, weighted by duration
    
    display(['Plate #' num2str(cur_file_num) ' : ']);
    display(['   IT N: ' num2str(all_IT_N)]);
    display(['   IT times (sec): ' num2str(mean(all_IT_t)) ' p/m ' num2str(std(all_IT_t)/length(all_IT_t))]);
    display(['   IT sum t (sec): ' num2str(all_IT_sum_t)]);
    display(['   IT fraction (time>MINTIME_IT) : ' num2str(all_IT_time_frac)]);
    display(['   IT fraction (dist) : ' num2str(all_IT_dist_frac)]);
    display(['   IT band: ' num2str(all_IT_av_T) '^oC \pm ' num2str(all_IT_std_T) ]);
    display('-------');
    IT_t = [IT_t; all_IT_t];
    IT_N = [IT_N; all_IT_N];
    IT_sum_t = [IT_sum_t; all_IT_sum_t];
    IT_time_fracs = [IT_time_fracs; all_IT_time_frac];
    IT_dist_fracs = [IT_dist_fracs; all_IT_dist_frac];
    IT_center_T = [IT_center_T; all_IT_av_T];
    IT_width_T = [IT_width_T; all_IT_std_T];
    
    figure(5000+k); 
    title(['Sample tracks, plate #' num2str(cur_file_num)]);
    next_tr_with_IT = 0;
    for j = 1:9
        subplot(3,3,j);
        next_tr_with_IT = next_tr_with_IT+1;
        while next_tr_with_IT < length(cur_tr) & cur_tr(next_tr_with_IT).IT_num == 0
        	next_tr_with_IT = next_tr_with_IT+1;
        end;

        if next_tr_with_IT < length(cur_tr) & cur_tr(next_tr_with_IT).IT_num > 0   % next_tr_with_IT not out of bounds   
            i1 = cur_tr(next_tr_with_IT).IT_indx(1,1); % start IT
            i2 = cur_tr(next_tr_with_IT).IT_indx(1,2); % end IT
            plot(cur_tr(next_tr_with_IT).x(i1:i2),cur_tr(next_tr_with_IT).y(i1:i2));
            axis('equal');
        end;
        
    end; % for j=1:9
end; % for p = 1:length(plate_nums)
runlengthanal4(tr,plate_nums); % for plots with all plates

if length(plate_nums)>1 % more than one plate
    display(['All plates: ']);
    display(['   IT total number: ' num2str(mean(IT_N)) ' p/m ' num2str(std(IT_N)/length(IT_N))]);
    display(['   IT times (sec): ' num2str(mean(IT_t)) ' p/m ' num2str(std(IT_t)/length(IT_t))]);
    display(['   IT sum t (sec): ' num2str(mean(IT_sum_t)) ' p/m ' num2str(std(IT_sum_t)/length(IT_sum_t))]);
    display(['   IT fraction (time>MINTIME_IT) : ' num2str(mean(IT_time_fracs)) ' p/m ' num2str(std(IT_time_fracs)/length(IT_time_fracs))]);
    display(['   IT fraction (dist) : '  num2str(mean(IT_dist_fracs)) ' p/m ' num2str(std(IT_dist_fracs)/length(IT_dist_fracs))]);
    display(['   IT band center: '  num2str(mean(IT_center_T)) ' p/m ' num2str(std(IT_center_T)/length(IT_center_T))]);
    display(['   IT band width: '  num2str(mean(IT_width_T)) ' p/m ' num2str(std(IT_width_T)/length(IT_width_T))]);
    display('-------');

    I = num2str(mean(Icr_by_plate));
    Ie = num2str(std(Icr_by_plate)/sqrt(length(Icr_by_plate)));
    N = num2str(mean(Ns_by_plate,1));
    Ne = num2str(std(Ns_by_plate,1)/sqrt(size(Ns_by_plate,1)));
    D = num2str(mean(Durs_by_plate,1));
    De = num2str(std(Durs_by_plate,1)/sqrt(size(Durs_by_plate,1)));
    V = num2str(mean(Vs_by_plate,1));
    Ve = num2str(std(Vs_by_plate,1)/sqrt(size(Vs_by_plate,1)));
  
    %%% Analyze all at once %%%    
    % [Icr,Ns,Durs, Vs] = runlengthanal4(tr,plate_nums);
    display(['I_c = ' I '+/-' Ie]);
    display(['Ns = ' N '+/-' Ne]);
    display(['Durs = ' D '+/-' De]);
    display(['Vs = ' V '+/-' Ve]);
    multiple_plates_meanT_over_time(tr);
end;

display(['Run time: ' num2str(round(etime(clock,t0)/6)/10) ' minutes.']);  % time in minutes

