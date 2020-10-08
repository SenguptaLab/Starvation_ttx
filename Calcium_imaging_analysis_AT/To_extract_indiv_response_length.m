%This runs on each blcct file

sample_rate = 0.5;
prompt = {'how long is the experiment in frames?'};
answer_length = inputdlg(prompt);
analysis_length = str2double(answer_length{1});
button = questdlg('Do you want to analyze on BLC_raw_delta_F_partial?','data sets','Yes','No','No');
answer = strcmp(button, 'Yes');
 if answer==0;
        BLC_raw_delta_F = BLC_raw_delta_F(1:analysis_length);
 elseif answer==1;
     BLC_raw_delta_F = BLC_raw_delta_F_partial(1:analysis_length);
 end
baseline_Ca_value = 10;
Image_Time = 1:analysis_length;
    for p = 1:analysis_length
         
        if  BLC_raw_delta_F(p) > baseline_Ca_value
            time_above_baseline(p) = 0.5;
            else
            time_above_baseline(p) = NaN;
         end
         
        
    end
    
   index = find(~isnan(time_above_baseline));
   T = find(diff(index)~=1)>=1;
   if T ==1
   idx = find(diff(index)~=1);
   A = [idx(1) diff(idx) numel(index)-idx(end)];
   else
       A = max(index)-min(index)+1;
   end
   each_response_duration = A'/2
   average_duration = mean(each_response_duration)
   total_response = sum(each_response_duration)
   clear all
   