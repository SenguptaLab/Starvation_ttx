function [] = make_plot_turns_in_time_bins(ds)

turns_per_min = ds.plates.turns_num_in_intrvl/5;
s = sum(turns_per_min,2); 
intrvls = size(turns_per_min,1);
centers = 2.5:5:(5*intrvls); % each row is a 5 minute intervals
figure;
clf;
bar(centers, s, 0.9, 'facecolor', 0.9*[1 1 1.07]); 
title('Sum numbers of turns from all plates');
xlabel('time (min)')
ylabel('Average number of turns per minute (all worms on plate)');

return;