clear; close all
dir_of_files = uigetdir();
addpath(dir_of_files);
D = dir([dir_of_files filesep '*.mat']);
GC = [];
baseline_ca_value = 10;
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
    %all_traces = vertcat(all_traces, BLC_raw_delta_F_partial(1:exp_length));
   
end

bins_response = [];
%bins = 1:20:301;
bins = 1:30:exp_length;
edge_end = (exp_length - 1)/30;
edges = [0, 1:1:edge_end edge_end];
%edges = [0, 1:1:10 10];
proportion = zeros(length(D), length(bins)-1);
for l= 1:length(bins)-1;
     bins_response = all_traces(:,bins(l):bins(l+1));
     for m = 1:length(D);
        proportion (m, l)= max(bins_response(m, :));
     end
end
proportion>=baseline_ca_value;
histgram_proportion = [];
for n= 1:length(bins)-1;
    histgram_proportion(1, n)= sum(ans(:,n))/length(D);
end
bar(histgram_proportion,1, 'edgecolor', [1 1 1], 'Linewidth', 0.5);
hold on
%xt = {'0'; '30'; '60'; '90'; '120'; '150'; '180'; '210'; '240'; '270'; '300'};
xt = {'0'; '20'; '40'; '60'; '80'; '100'; '120'; '140'; '160'; '180'; '200'; '220'; '240'; '260'; '280'; '300'};
set(gca, 'xticklabel', xt);
ylim([0 1]);
      
        