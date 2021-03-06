function [] = make_plot_run_distributions2(ds)

global ANG_FOR_RW

[N1,centers1,N3,centers3] = bin_runs_up_down4(ds); % determine centers for all plates

dn_by_plt = []; % each column - a plate, each row - a bin of down-runs durations
up_by_plt = []; % each column - a plate, each row - a bin of up-runs durations
pd = 0; % number  of plates with down-runs
pu = 0; % number  of plates with up-runs
for k = 1:length(ds.runs.time) % count ITs in each bin for each plate with ITs
    ang = ds.runs.angle{k};
    ang(ang>pi) = 2*pi-ang(ang>pi); % there shouldn't be any angels >pi but..   
    ang = abs(ang); % in terms of up/down gradient alpha is equivalent to -alpha 
    indx_dn = find(ang>(pi-ANG_FOR_RW) & ang<=pi);
    indx_up = find(ang<ANG_FOR_RW & ang>=0);
    runs_by_plt = ds.runs.time{k};
    if ~isempty(indx_dn) % only count plates with runs-down
        pd = pd+1;       % one more plate with runs-down
        N = hist(runs_by_plt(indx_dn),centers1); % bin runs-down in plate
        dn_by_plt(:,pd) = N'; % add binned runs to "by-plate" matrix
    end;
    if ~isempty(indx_up) % only count plates with ITs
        pu = pu+1;       % one more plate with runs-up
        N = hist(runs_by_plt(indx_up),centers3); % bin runs-up in plate
        up_by_plt(:,pu) = N'; % add binned runs to "by-plate" matrix
    end;
end;

N1 = mean(dn_by_plt, 2);
eN1 = std(dn_by_plt,1,2)/size(dn_by_plt,2);
NN1 = 100*N1/sum(N1); % convert to percentage
eNN1 = 100*eN1/sum(N1);
N3 = mean(up_by_plt, 2);
eN3 = std(up_by_plt,1,2)/size(up_by_plt,2);
NN3 = 100*N3/sum(N3); % convert to percentage
eNN3 = 100*eN3/sum(N3);

indx1 = find(N1>0);
if ~isempty(indx1) 
    [A1, lambda1] = my_exp_fit4(centers1', N1);
    [NA1, Nlambda1] = my_exp_fit4(centers1', NN1);
else
    A1 = 0;
    lambda1 = 0; 
    NA1 = 0;
    Nlambda1 = 0; 
end;

indx3 = find(N3>0);
if ~isempty(indx3) 
    [A3, lambda3] = my_exp_fit4(centers3', N3);
    [NA3, Nlambda3] = my_exp_fit4(centers3', NN3);
else
    A3 = 0;
    lambda3 = 0; 
    NA3 = 0;
    Nlambda3 = 0; 
end;

% Non-normalized runs UP/DOWN histogram
fh = figure;
clf;
set(fh,'position',[621   400   756   567]);
bar(centers1, N1, 0.90, 'c'); % down
hold on;
for k = 1:length(N1)
    line(centers1(k)+[0 0], N1(k)+eN1(k)*[-1 1], 'color', 0.4*[1 1 1]);
end;
bar(centers3, N3, 0.78, 'r'); % up
for k = 1:length(N3)
    line(centers3(k)+[0 0], N3(k)+eN3(k)*[-1 1], 'color', 0.4*[1 1 1]);
end;
plot(centers1(indx1), A1*exp(-1*centers1(indx1)/lambda1), '--','linewidth',2,'color',[0 0.4 1]);
plot(centers3(indx3), A3*exp(-1*centers3(indx3)/lambda3), '--','linewidth',2,'color',[1 0.4 0]);
hold off;
title(['Non-Normalized run durations (averaged over plates), \tau_{dn} = ' num2str(lambda1) ' \tau_{up} = ' num2str(lambda3) ]);
ch = get(fh,'children');
set(ch,'yscale','log');
xlabel(['Mean duration of binned runs (sec)']);
ylabel(['Number of runs']);

%%% also calculate intersection point of two fits:
%tstar = (lambda1*lambda3)/(lambda1-lambda3) * log(A3/A1);
%ystar = A1*exp(-1*tstar/lambda1); 
%hold on;
%plot(tstar,ystar,'g*');
%hold off;
%text(tstar,ystar*3,['t^* = ' num2str(tstar)]);


% Normalized runs UP/DOWN histogram
fh = figure;
clf; 
set(fh,'position',[121   100   672   504]);
bar(centers1, NN1, 0.90, 'c');
hold on; 
for k = 1:length(NN1)
    line(centers1(k)+[0 0], NN1(k)+eNN1(k)*[-1 1], 'color', 0.4*[1 1 1]);
end;
bar(centers1, NN3,0.78,'r');
for k = 1:length(NN3)
    line(centers3(k)+[0 0], NN3(k)+eNN3(k)*[-1 1], 'color', 0.4*[1 1 1]);
end;
plot(centers1(indx1), NA1*exp(-1*centers1(indx1)/Nlambda1), '--','linewidth',2,'color',[0 0.4 1]);
plot(centers3(indx3), NA3*exp(-1*centers3(indx3)/Nlambda3), '--','linewidth',2,'color',[1 0.4 0]);
title(['Normalized run durations (averaged over plates), \tau_{dn} = ' num2str(Nlambda1) ' \tau_{up} = ' num2str(Nlambda3)]);
hold off;
ch = get(fh,'children');
set(ch,'yscale','log');
xlabel(['Mean duration of binned runs (sec)']);
ylabel(['Percentage of runs (from total number of runs)']);

return;