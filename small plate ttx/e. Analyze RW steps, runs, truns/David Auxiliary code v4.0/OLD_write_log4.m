function [Icr,Ns,Durs,Vs]=write_log4(tr, res, plate_nums, plates_indx)

global FPS; 

dur_dn = res.dur_dn;
vel_dn = res.vel_dn;
crv_dn = res.crv_dn;
vel_crv_dn = res.vel_crv_dn;

dur_vr = res.dur_vr;
vel_vr = res.vel_vr;
crv_vr = res.crv_vr;
vel_crv_vr = res.vel_crv_vr;

dur_up = res.dur_up;
crv_up = res.crv_up;
vel_up = res.vel_up;
vel_crv_up = res.vel_crv_up; 

%%% Make (and display) log file %%%
dir_name = tr(1).path; 
tmp = strfind(dir_name,'\');
if length(plate_nums)==1 
    pref = ['LOG_p' num2str(plate_nums(1)) '_'];
    disp(['Logging plate #' num2str(plate_nums(1))]);
else
    pref = ['LOG_all_'];
	disp(['Logging plates ' num2str(plate_nums)]);
end;
sv_nm = [pref dir_name((tmp(end)+1):end) '_'];
tmp = datestr(now,31);
tmp = ['D' tmp(1:end-3)]; % get rid of seconds
tmp = regexprep(tmp,{':'},'_');
tmp = regexprep(tmp,{' '},'_H');
sv_nm = [sv_nm tmp '.txt'];

lh = fopen([dir_name '\' sv_nm],'w');

%%% Total durations (run_num x run_duration) in each direction
str = 'sum total duration of runs (sec): ';
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
str = ['dn=' num2str(sum(dur_dn)) ' ; vr=' num2str(sum(dur_vr)/2)];
str = [str ' ; up=' num2str(sum(dur_up))];
str = [str ' (' num2str(100*sum(dur_dn)/sum(dur_up)-100) '%)'];
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
fprintf(lh,'\r\n');

%%% Total duration of runs down without long-runs tail 
if ~isempty(dur_up)
    mx_up = 0.8*max(dur_up);
else
    mx_up = 0.6*max(dur_dn);
end;
indx1 = find(dur_dn<mx_up); % all head
indx2 = find(dur_dn>=mx_up); % all tail 
frac2keep = max(0,length(dur_up)-length(indx2))/length(indx1); % required size of head: (#up-#tail)/#head
AllDown = sum(dur_dn);
AllUp = sum(dur_up);
AllHead = sum(dur_dn(indx1));
AllTail = sum(dur_dn(indx2));
ExcessHead = round(AllHead*(1-frac2keep));
str = 'Runs-down durations: All, Tail, Excess head (sec)';
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
str = [num2str(AllDown) ' ' num2str(AllTail) ' ' num2str(ExcessHead)];
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
str = 'Bias in total migration: All, No tail, No Excess head (sec)';
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
ci1 = (AllDown-AllUp)/(AllDown+AllUp);
ci2 = (AllDown-AllTail-AllUp)/(AllDown-AllTail+AllUp);
ci3 = (AllDown-ExcessHead-AllUp)/(AllDown-ExcessHead+AllUp);
str = [num2str(ci1) ' ' num2str(ci2) ' ' num2str(ci3)];
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');

%%% number of runs in each direction
Ns = [length(dur_dn), length(dur_vr)/2, length(dur_up)];
str = '# of runs dn/vr/up:';
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
str = ['dn=' num2str(Ns(1)) ' ; vr=' num2str(Ns(2)/2)];
str = [str ' ; up=' num2str(Ns(3))];
str = [str ' (' num2str(100*Ns(1)/Ns(3)-100) '%)'];
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
fprintf(lh,'\r\n');

%%% number of curved runs in each direction
N_C_s = [length(crv_dn), length(crv_vr)/2, length(crv_up)];
str = '# of CURVED runs:';
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
str = ['dn=' num2str(N_C_s(1)) ' ; vr=' num2str(N_C_s(2)/2)];
str = [str ' ; up=' num2str(N_C_s(3))];
str = [str ' (' num2str(100*N_C_s(1)/N_C_s(2)-100) '%)'];
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
fprintf(lh,'\r\n');

% Duration of runs in each direction
Durs = [mean(dur_dn), mean(dur_vr) ,mean(dur_up)];
str = '<duration> of runs dn/vr/up (sec):';
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
str = ['dn=' num2str(Durs(1)) ' ; vr=' num2str(Durs(2))];
str = [str ' ; up=' num2str(Durs(3))];
str = [str ' (' num2str(100*Durs(1)/Durs(3)-100) '%)'];
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
fprintf(lh,'\r\n');

% Tau of runs in each direction (from fit to exponential distributions)
str = '\tau (from normalized hists) of runs dn/vr/up (sec):';
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
str = ['dn=' num2str(res.Ntau_dn) ' ; vr=' num2str(res.Ntau_vr)];
str = [str ' ; up=' num2str(res.Ntau_up)];
str = [str ' (' num2str(100*res.Ntau_dn/res.Ntau_up-100) '%)'];
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
fprintf(lh,'\r\n');

% Duration of curved runs in each direction
Dur_C_s = [mean(crv_dn), mean(crv_vr) ,mean(crv_up)];
str = '<duration> of CURVED runs (sec):';
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
str = ['dn=' num2str(Dur_C_s(1)) ' ; vr=' num2str(Dur_C_s(2))];
str = [str ' ; up=' num2str(Dur_C_s(3))];
str = [str ' (' num2str(100*Dur_C_s(1)/Dur_C_s(3)-100) '%)'];
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
fprintf(lh,'\r\n');

% velocities of runs in each direction
Vs = [mean(vel_dn), mean(vel_vr), mean(vel_up)];
str = '<velocity> during runs (pix/sec):';
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
str = ['dn=' num2str(Vs(1)) ' ; vr=' num2str(Vs(2))];
str = [str ' ; up=' num2str(Vs(3))];
str = [str ' (' num2str(100*Vs(1)/Vs(3)-100) '%)'];
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
fprintf(lh,'\r\n');

% velocities of curved runs in each direction
V_C_s = [mean(vel_crv_up), mean(vel_crv_vr), mean(vel_crv_dn)];
str = '<velocity> during CURVED runs (pix/sec):';
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
str = ['dn=' num2str(V_C_s(1)) ' ; vr=' num2str(V_C_s(2))];
str = [str ' ; up=' num2str(V_C_s(3))];
str = [str ' (' num2str(100*V_C_s(1)/V_C_s(3)-100) '%)'];
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
fprintf(lh,'\r\n');

Icr = (sum(dur_dn)-sum(dur_up))/(sum(dur_dn)+sum(dur_up));
Icr_av  = (mean(dur_dn)-mean(dur_up))/(mean(dur_dn)+mean(dur_up));
I_tau = (res.tau_dn-res.tau_up)/(res.tau_dn+res.tau_up);
str = 'I_cr:';
str = [str 'tot=' num2str(Icr) ' ; av=' num2str(Icr_av)];
str = [str ' ; \tau=' num2str(I_tau)];
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
fprintf(lh,'\r\n');

wormsec=0;
for i_indx = 1:length(plates_indx) % for each track of each worm (from allowed plates)
    i = plates_indx(i_indx);
    wormsec=wormsec+length(tr(i).t);
end
str = ['total worm hours tracked = ' num2str(wormsec/3600)];
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
fprintf(lh,'\r\n');

allRuns = 0;
allITs = 0;
for k = 1:length(tr)
    allRuns = allRuns+length(tr(k).run_angle);
    allITs = allITs+length(tr(k).IT_angle);
end;
str = ['percentage of IT: ' num2str(allITs/allRuns*100) '%'];  % time in minutes
fwrite(lh, str, 'char');
fprintf(lh,'\r\n');
fprintf(lh,'\r\n');

% IT info


fclose(lh); 
str = textread([dir_name '\' sv_nm], '%s', 'delimiter', '\n');
disp(str);
fclose('all');

return;