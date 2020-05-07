% This script will plot a heatmap from any group of blccrt or GCaMP_Analysis
% or FRET_Analysis files. It will prompt you for how long your plot should be, in frames.
% FRAMES, NOT SECONDS. Enter the number of frames you want plotted and a
% heat map of that length will be plotted. It does some minor smoothing for
% prettiness.
clear; close all
dir_of_files = uigetdir();
addpath(dir_of_files);
D = dir([dir_of_files filesep '*.mat']);
GC = [];

prompt = {'How long is your experiment(in frames)?:'};
dlg_title = 'Analysis setting';
num_lines = 1;
def = {'241'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
exp_length = str2double(answer{1});
button = questdlg('Enter the upper threshold of heatmap:','Threshold','30','50','other', '30');
      switch button
          case '30';
              Upper_threshold = 30;
          case '50'
              Upper_threshold = 50;
          case 'other'
              prompt = {'Enter the upper threshold of heatmap:'};
              answer_threshold = inputdlg(prompt);
              Upper_threshold = str2double(answer_threshold{1}); 
      end
all_traces = [];
%all_temp = [];
for i = 1:length(D);
    BLC_file_name = D(i).name(1:end-4);
    load([dir_of_files filesep BLC_file_name]);
   
    all_traces = vertcat(all_traces, BLC_raw_delta_F(1:exp_length));
    %all_temp = vertcat(all_temp, Temperature(1:exp_length));
   
end

colormap('hot')
h = imagesc(all_traces,[0 Upper_threshold]);
xax = get(gca,'XTick');
rt = xax / 2;
set(gca,'XTicklabel',rt)
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
cb = colorbar;
date_proc = datevec(now);

saveas(gcf,[dir_of_files filesep 'heatmap.fig']);
saveas(gcf,[dir_of_files filesep 'heatmap.jpg'],'jpg')
print(gcf,[dir_of_files filesep 'heatmap.eps'],'-depsc2')

pause_time = 3;
pause(pause_time);

%TO sort
maxAll_traces = max(all_traces, [], 2);
[max_signal, index] = sort(maxAll_traces, 'descend');
sort_traces = all_traces(index, :);
colormap('jet')
h = imagesc(sort_traces,[0 Upper_threshold]);
xax = get(gca,'XTick');
rt = xax / 2;
set(gca,'XTicklabel',rt)
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
cb = colorbar;
date_proc = datevec(now);

saveas(gcf,[dir_of_files filesep 'heatmap_sort.fig']);
saveas(gcf,[dir_of_files filesep 'heatmap_sort.jpg'],'jpg')
print(gcf,[dir_of_files filesep 'heatmap_sort.eps'],'-depsc2')

%TO sort in the order of timing of first response
baseline_ca_value = 10;
sz = [length(D) 1];
non_responder = exp_length + 1;
index_timing = ones(sz)*non_responder;
for j= 1:length(D);
    if max(all_traces(j, 5:end))<10,
       index_timing(j, 1)= non_responder;
    else    
    index_timing(j, 1) = find(all_traces(j, 5:end)>baseline_ca_value, 1, 'first');
    end
end
[timing_signal, index2] = sort(index_timing, 'ascend');
sort_traces_timing = all_traces(index2, :);
colormap('jet')
h = imagesc(sort_traces_timing,[0 Upper_threshold]);
xax = get(gca,'XTick');
rt = xax / 2;
set(gca,'XTicklabel',rt)
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
cb = colorbar;
date_proc = datevec(now);

saveas(gcf,[dir_of_files filesep 'heatmap_sort_timing.fig']);
saveas(gcf,[dir_of_files filesep 'heatmap_sort_timing.jpg'],'jpg')
print(gcf,[dir_of_files filesep 'heatmap_sort_timing.eps'],'-depsc2')

pause(pause_time);

