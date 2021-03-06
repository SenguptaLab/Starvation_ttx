function [Icr,Ns,Durs,Vs]=runlengthanal4(tr,varargin)

% INPUT:
%   tr - a structure pepared by flagturns3.m containing tracks and turns
%        data (see flagturn3.m)
%   [plates] - an optional list of plate numbers to analyze
% PARAMETERS:
%   FPS - frames per second (global) 
% OUTPUT:
%   cents - 
%   mdur - 
%   sdur - 


global TIME_BIN
global ANG_BIN_1
global MAX_RUN_DURATION_FOR_RW
global MIN_RUN_DURATION_FOR_RW
global ANG_FOR_RW

warning off MATLAB:divideByZero
PLOTHISTS = 0;

%%% Find indices of cells in tr2 that belong to plates %%%
%%% that were chosen for analysis                      %%%
%%% Empty plate list means analyze ALL plates          %%%
if nargin==2       % second argument has plate no. list
    plates_indx = [];
    plate_nums = [];
    for k=1:length(tr)
        cur_id = tr(k).plate_id; 
        flag = ismember(cur_id,varargin{1}) || isempty(varargin{1});
        if  flag % find index in tr2 of tracks from listed plates (or all if list is empty)
            plates_indx = [plates_indx, k];
            plate_nums = [plate_nums, cur_id];
        end;
    end;
else  % all of the tracks in tr2 are to be analyzed
    plates_indx = 1:length(tr); 
end;
plate_nums = unique(plate_nums);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Analyze RW runs, turns, angles etc. 
ang = []; % average andle of each run
dur = []; % in seconds - actual run time for each run
allx = []; % total pixel moved in X direction for each run
ally = []; % total pixel moved in Y direction for each run
len = []; % distance in pixels (along contour of each run)
curvd = [];
for i_indx = 1:length(plates_indx) % for each track of each worm (from allowed plates)
    i = plates_indx(i_indx); 
    t = tr(i);
    indx = find((t.run_time)<MAX_RUN_DURATION_FOR_RW  &  (t.run_time)>MIN_RUN_DURATION_FOR_RW);
  
    allx = [allx; t.run_dx(indx)];         % end-to-end x distance for each run
    ally = [ally; t.run_dy(indx)];         % end-to-end y distance for each run
    dur  = [dur;  t.run_time(indx)];       % no. of frames for each run
    len  = [len;  t.run_dist(indx)];       % distance in pixels (along contour of each run)
    ang  = [ang;  t.run_angle(indx)];      % average angle in radians for each run
    curvd = [curvd; t.curved(indx)]; 
end;
vel = sqrt(allx.^2 + ally.^2)./dur;

% group symmetric angles so all are in range [0,pi]
ang(ang>pi)=2*pi-ang(ang>pi); % there shouldn't be any angels >pi but..   
ang=abs(ang); % in terms of up/down gradient alpha is equivalent to -alpha 
disp([num2str(length(dur)) ' runs found']);

timebin = TIME_BIN;
timeedges=[timebin:timebin:max(dur)]; % time-bin edges
timecents=[3*timebin/2:timebin:timeedges(end)]; % time-bin centers

bin1 = ANG_BIN_1;
edges1=[0:bin1:pi];
cents1=edges1(1)+bin1/2:bin1:edges1(end);
[n,inds1]=histc(ang,edges1);

u=1:length(edges1);
for i=1:length(u)-1;
    choose=find(inds1==u(i));
    angledurpoints{i}=dur(choose);
    mlen(i)=mean(len(choose));
    mdur(i)=mean(dur(choose));
    slen(i)=std(len(choose))/sqrt(length(choose));
    sdur(i)=std(dur(choose))/sqrt(length(choose));
    if PLOTHISTS
        n=histc(dur(choose),timeedges);
        n=n(1:end-1);
        [p,s]=lsq('modelexp',[100 10],timecents(n>0)',n(n>0),sqrt(n(n>0)),1);
        set(gca,'yscale','log','xlim',[min(timeedges),max(timeedges)]);
        fdur(i)=p(2);
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time 
indx_dn = find(ang>(pi-ANG_FOR_RW) & ang<=pi);
indx_vr = find(ang<(pi/2+ANG_FOR_RW) & ang>pi/2-ANG_FOR_RW);
indx_up = find(ang<ANG_FOR_RW & ang>=0);
indx_curvd = find(curvd==1); 

x_dn   = allx(indx_dn);
y_dn   = ally(indx_dn);
dur_dn = dur(indx_dn);
vel_dn = vel(indx_dn);
crv_dn = dur(intersect(indx_dn,indx_curvd));
vel_crv_dn = vel(intersect(indx_dn,indx_curvd));
res.dur_dn = dur_dn;
res.vel_dn = vel_dn;
res.crv_dn = crv_dn;
res.vel_crv_dn = vel_crv_dn;

x_vr = allx(indx_vr);
y_vr = ally(indx_vr);
dur_vr = dur(indx_vr);
vel_vr = vel(indx_vr);
crv_vr = dur(intersect(indx_vr,indx_curvd));
vel_crv_vr = vel(intersect(indx_vr,indx_curvd));
res.dur_vr = dur_vr;
res.vel_vr = vel_vr;
res.crv_vr = crv_vr;
res.vel_crv_vr = vel_crv_vr;

x_up   = allx(indx_up);
y_up   = ally(indx_up);
dur_up = dur(indx_up);
vel_up = vel(indx_up);
crv_up = dur(intersect(indx_up,indx_curvd));
vel_crv_up = vel(intersect(indx_up,indx_curvd));
res.dur_up = dur_up;
res.vel_up = vel_up;
res.crv_up = crv_up;
res.vel_crv_up = vel_crv_up;

%%% Calculate Ns, centers and taus for histograms %%%
num_of_bins = round(max(dur_dn)/8); % make ~8 seconds bins in the histograms
[N1, centers1] = hist(dur_dn,num_of_bins);
NN1 = 100*N1/sum(N1); % convert to percentage
if ~isempty(find(N1>0)) 
    [A1, lambda1] = my_exp_fit4(centers1, N1);
    [NA1, Nlambda1] = my_exp_fit4(centers1, NN1);
    res.tau_dn = lambda1; 
    res.Ntau_dn = lambda1; 
else
    res.tau_dn = 0; 
    res.Ntau_dn = 0; 
end;

[N2, centers2] = hist(dur_vr,centers1); % use same bins as down
NN2 = 100*N2/sum(N2); % convert to percentage
if ~isempty(find(N2>0)) 
    [A2, lambda2] = my_exp_fit4(centers2, N2);
    [NA2, Nlambda2] = my_exp_fit4(centers2, NN2);
    res.tau_vr = lambda2; 
    res.Ntau_vr = Nlambda2; 
else
    res.tau_vr = 0; 
    res.Ntau_vr = 0; 
end;

[N3, centers3] = hist(dur_up,centers1); % use same bins as down
NN3 = 100*N3/sum(N3); % convert to percentage
if ~isempty(find(N3>0)) 
    [A3, lambda3] = my_exp_fit4(centers3, N3);
    [NA3, Nlambda3] = my_exp_fit4(centers3, NN3);
    res.tau_up = lambda3;
    res.Ntau_up = Nlambda3;
else
    res.tau_up = 0; 
    res.Ntau_up = 0; 
end;

[Icr,Ns,Durs,Vs]=write_log4(tr,res,plate_nums,plates_indx);
run_err = [];
run_dist = [];
for k = 1:length(tr)
    run_err = [run_err; tr(k).run_err];
    run_dist = [run_dist; tr(k).run_dist];
end;
indx = find(run_dist>0); % should be all but just in case
disp(['Mean run relative error: ' num2str(mean(run_err(indx)./run_dist(indx)))]);

if length(plate_nums)<2 % do not make figs with one plate.. 
    return; 
end;

%%% Make figures %%%
SHAPE='o';  % shape in plot
figure(101); % run duration vs. angle
clf;
hold on;
%centsfine=[0:.01:pi];
[pf,sf]=lsq('modelcos',[1 1],cents1,mdur,sdur,0);
[pf2,sf2]=lsq('modelsgncos',[1 1],cents1,mdur,sdur,0);  % sign of cos fit...
[pl,sl]=lsq('modellinear',[1 1],cents1,mdur,sdur,0);
[pc,sc]=lsq('modelconst',[30],cents1,mdur,sdur,0);
h1=plot(cents1*180/pi,mdur,[SHAPE 'k-']); 
for i=1:length(cents1)
    h1=line(180/pi*[cents1(i),cents1(i)],[mdur(i)-sdur(i),mdur(i)+sdur(i)]); 
    set(h1,'color','k');
end;
title([': run duration vs. angle, P_{const}=' num2str(sc.p_fit) ', P_{cos}= ' num2str(sf.p_fit) ...
        ', P_{sgn(cos)}= ' num2str(sf2.p_fit) ', P_{lin}= ' num2str(sl.p_fit)]);
ylabel('mean run duration (seconds)');
xlabel(['angle center']);
set(gca,'ylim',[10,40]);

figure(102); % runs from common origin in pixels
clf;
plot(allx,ally,'k.'); hold on;
set(gca,'dataaspect',[1 1 1]);
h1=plot(0,0,'.');
set(h1,'color','white');
title('runs from common origin in pixels');

figure(103); % runs from common origin in seconds
plot(dur_up.*x_up./sqrt(x_up.^2+y_up.^2),dur_up.*y_up./sqrt(x_up.^2+y_up.^2),'r.'); 
hold on;
plot(dur_dn.*x_dn./sqrt(x_dn.^2+y_dn.^2),dur_dn.*y_dn./sqrt(x_dn.^2+y_dn.^2),'b.'); 
plot(dur_vr.*x_vr./sqrt(x_vr.^2+y_vr.^2),dur_vr.*y_vr./sqrt(x_vr.^2+y_vr.^2),'k.'); 
h1=plot(0,0,'.');
set(h1,'color','white');
title(['runs from common origin in seconds']);
xlabel('time (s)');
ylabel('time (s)');
axis 'equal';

fh1 = figure(201); % Normalized histogram of durations of runs DOWN the gradient
indx1 = find(NN1>0);
bar(centers1, NN1);
hold on;
plot(centers1(indx1), NA1*exp(-1*centers1(indx1)/Nlambda1), 'r--','linewidth',2);
hold off;
title(['durations of runs ~DOWN the gradient, \tau = ' num2str(Nlambda1)]);
ch1 = get(fh1,'children');
set(ch1,'yscale','log');

fh2 = figure(202); % Normalized histogram of durations of runs VERTICAL to the gradient
indx2 = find(NN2>0);
bar(centers2, NN2);
hold on;
plot(centers2(indx2), NA2*exp(-1*centers2(indx2)/Nlambda2), 'r--','linewidth',2);
hold off;
title(['durations of runs ~VR to the gradient, \tau = ' num2str(Nlambda2)]);
ch2 = get(fh2,'children');
set(ch2,'yscale','log');

fh3 = figure(203); % Normalized histogram of durations of runs UP the gradient
bar(centers3, NN3);
hold on;
indx3 = find(NN3>0);
plot(centers3(indx3), NA3*exp(-1*centers3(indx3)/Nlambda3), 'r--','linewidth',2);
hold off;
title(['durations of runs ~UP the gradient, \tau = ' num2str(Nlambda3)]);
ch3 = get(fh3,'children');
set(ch3,'yscale','log');

figure(204); % NON-Normalized histograms of durations of runs DOWN/UP 
close(204);
fh4 = figure(204);
set(fh4,'position',[121   100   840   630]);
bar(centers1, N1, 'c'); % down
hold on;
bar(centers3, N3, 'r'); % up
plot(centers1(indx1), A1*exp(-1*centers1(indx1)/lambda1), '--','linewidth',2,'color',[0 0.4 1]);
plot(centers3(indx3), A3*exp(-1*centers3(indx3)/lambda3), '--','linewidth',2,'color',[1 0.4 0]);
hold off;
title(['Non-Norm run durations, \tau_{dn} = ' num2str(lambda1) ' \tau_{up} = ' num2str(lambda3) ]);
ch4 = get(fh4,'children');
set(ch4,'yscale','log');
%%% also calculate intersection point of two fits:
tstar = (lambda1*lambda3)/(lambda1-lambda3) * log(A3/A1);
ystar = A1*exp(-1*tstar/lambda1); 
hold on;
plot(tstar,ystar,'g*');
hold off;
text(tstar,ystar*3,['t^* = ' num2str(tstar)]);

fh5 = figure(205); % Normalized histogram of durations of runs DOWN/UP
clf; 
bar(centers1, 100*N1/sum(N1),0.9,'c');
hold on; 
bar(centers1, 100*N3/sum(N3),0.67,'r');
plot(centers1(indx1), NA1*exp(-1*centers1(indx1)/Nlambda1), '--','linewidth',2,'color',[0 0.4 1]);
plot(centers3(indx3), NA3*exp(-1*centers3(indx3)/Nlambda3), '--','linewidth',2,'color',[1 0.4 0]);
title(['Norm run durations, \tau_{dn} = ' num2str(Nlambda1) ' \tau_{up} = ' num2str(Nlambda3)]);
hold off;
ch5 = get(fh5,'children');
set(ch5,'yscale','log');

figure(301); % histogram of durations of CURVED runs DOWN the gradient
[N4, centers4] = hist(crv_dn);
%hist(dur_up);
bar(centers4, 100*N4/sum(N4));
title('durations of CURVED runs ~DOWN the gradient');

figure(302); % histogram of durations of CURVED runs VERTICAL to the gradient
[N5, centers5] = hist(crv_vr);
%hist(dur_dn);
bar(centers5, 100*N5/sum(N5));
title('durations of CURVED runs ~VERTICAL to the gradient');

figure(303); % histogram of durations of CURVED runs UP the gradient
[N6, centers6] = hist(crv_up);
%hist(dur_vr);
bar(centers6, 100*N6/sum(N6));
title('durations of CURVED runs ~UP the gradient');

figure(401); % fraction contribution of each bin to cryophilic drive 
bin_limits = [0 (centers1(2:end)+centers1(1:(end-1)))/2 999]; 
for k = 1:(length(bin_limits)-1)
    runs_in_bin_k = find(dur>=bin_limits(k) & dur<bin_limits(k+1));
    dn_bin_k = intersect(runs_in_bin_k,indx_dn);
    up_bin_k = intersect(runs_in_bin_k,indx_up);
    delta_by_bins(k) = sum(dur(dn_bin_k))-sum(dur(up_bin_k));
    runs_up_present(k) = sum(dur(up_bin_k))>0;
end;
plot(centers1,delta_by_bins/sum(delta_by_bins),'x:','markersize',8);
hold on;
plot(centers1(runs_up_present),delta_by_bins(runs_up_present)/sum(delta_by_bins),'or','markersize',8);
hold off;
title(['fraction contribution of each bin to cryophilic drive']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

