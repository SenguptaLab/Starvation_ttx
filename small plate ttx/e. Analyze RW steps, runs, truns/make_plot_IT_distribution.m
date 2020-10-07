function [] = make_plot_IT_distribution(ds)

BIN_WIDTH = 8; % seconds

ITs = [];
for k = 1:length(ds.runs.time)    
    ITs = [ITs; ds.ITs.time{k}];
end;

num_of_bins = round(max(ITs)/BIN_WIDTH); % make ~BIN_WIDTH seconds bins in the histograms
[N, centers] = hist(ITs,num_of_bins);
NN = 100*N/sum(N); % convert to percentage
indx = find(N>0);
if ~isempty(indx) 
    [A, lambda] = my_exp_fit4(centers, N);
    [NA, Nlambda] = my_exp_fit4(centers, NN);
else
    A = 0;
    lambda = 0; 
    NA = 0;
    Nlambda = 0; 
end;

% Non-normalized runs UP/DOWN histogram
fh = figure;
clf;
set(fh,'position',[621   400   756   567]);
bar(centers, N, 0.9, 'facecolor', 0.9*[1 1 1.07]); 
hold on;
plot(centers(indx), A*exp(-1*centers(indx)/lambda), '--','linewidth',2,'color',0.4*[1 1 1]);
hold off;
title(['Non-Normalized run durations, \tau_{IT} = ' num2str(lambda) ]);
ch = get(fh,'children');
set(ch,'yscale','log');
title('All ITs from all plates');
xlabel(['Mean duration of binned ITs (sec)']);
ylabel(['Number of runs']);

% Normalized runs UP/DOWN histogram
fh = figure;
clf; 
set(fh,'position',[121   100   672   504]);
bar(centers, NN, 0.9, 'facecolor', 0.9*[1 1 1.07]);
hold on; 
plot(centers(indx), NA*exp(-1*centers(indx)/Nlambda), '--','linewidth',2,'color',0.4*[1 1 1]);
title(['Normalized run durations, \tau_{NIT} = ' num2str(Nlambda) ]);
hold off;
ch = get(fh,'children');
set(ch,'yscale','log');
title('All ITs from all plates');
xlabel(['Mean duration of binned ITs (sec)']);
ylabel(['Percentage of ITs (from total number of runs)']);

return;
