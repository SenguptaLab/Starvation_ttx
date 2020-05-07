clear all 

c = clock;
path = ['C:\ThermoAssays\' num2str(c(2)) '.' num2str(c(1)) filesep];
fn = 1;
[fn,path] = uigetfile([path '*.mat'],'Select mat file; or click "cancel"');
load([path fn]);

baseline_partial=mean(Signal(181:190));
delta_F_partial=((Signal(181:481)-baseline_partial)/baseline_partial)*100;
Temperature_partial = Temperature(181:481);
Image_Time_partial =Image_Time(181:481);

figure(2)
[AX,H1,H2]=plotyy(Image_Time_partial,delta_F_partial, Image_Time_partial,Temperature_partial,'plot');
hold on
%[AX,H1,H2]=plotyy(Image_Time,FRET1, Image_Time,Temperature,'plot');
set(H1,'LineStyle','-','LineWidth',2,'color', 'b');
set(H2,'LineStyle','--','LineWidth',2,'color', 'g');
%xlim([181 481]);
%ylim(axes_Delta_R,[-50 150])
ylabel('\deltaF/F (%)');
xlabel('Time (sec) ');
ylabel(AX(2),'Temperature \circC')

saveas(gcf,[path 'Delta_F_Figure_partial_' ROI_name '.fig'],'fig');
saveas(gcf,[path 'Delta_F_Figure_partial_' ROI_name '.jpg'],'jpg');
Mat_file_name = ['GCAMP_Analysis_' ROI_name '.mat'];
save ([path 'partial_' Mat_file_name],'baseline','baseline_partial','BG_position','Signal','delta_F','delta_F_partial', 'Mat_file_name','h','i','Image_Time','Image_Time_partial', 'mypath','ROI','Norm_Signal','R_Min_Col','R_Neuron_BG_X1','R_Neuron_BG_X2','R_Neuron_BG_Y1','R_Neuron_BG_Y2','R_Neuron_X1','R_Neuron_X2', 'R_Neuron_Y1','R_Neuron_Y2','Right_Neuron_BG','Temperat','Temperature','Temperature_partial','ROI_name')
        %clear AX BG_position button ext figure h H1 H2 Image1 Images IX
        %clear L_Neuron_loc R_Neuron_Im R_Neuron_loc IX L_IX L_Neuron_BG
        %clear L_Neuron_Pixels R_Neuron_BG R_Neuron_Pixels T_fileR_Min_Col M 
        %clear AX Delta_R_Figure H1 H2 M Neuron_Pixel_Size axes_Delta_R
        clf
       
        if  exist([path filesep 'BLC'], 'dir') ==0;
             mkdir(path,'BLC');
        else 
             disp('BLC file added to the folder already exists');
        end
         set(0,'DefaultFigurePosition',[100 100 1000 1000]);
        fig(3)  = plot(delta_F_partial);
        disp(fn);
        title(fn,'Interpreter','none');
        ax =gca;
        axes(ax);
        [blx, bly] = getline(ax);
        [ublx, a, b] = unique(blx);
        ubly = bly(a);

        bsl = interp1(ublx,ubly,1:length(delta_F_partial));
        bsl(isnan(bsl)) = 0;

        close
        
        BLC_raw_delta_F_partial = (delta_F_partial + abs(min(bly))) - (bsl +  abs(min(bly)));

        figure;
        hold on
        plot(BLC_raw_delta_F_partial,'r');
        title('Corrected signal')
        pause(2)
        close

        BLC_file = (['blcrrt_partial_', ROI_name]);
        save([path, 'BLC/', BLC_file]);
close all 
clear all