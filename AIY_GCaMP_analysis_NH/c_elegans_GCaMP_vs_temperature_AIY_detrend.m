%This script incorporates deltaF/F calculations from Jenna Sternberg and
%members of the Wyart lab, and temperature measurements from Harry Bell in
%the Sengupta lab.

%Generates deltaF/F traces from multiple ROI in one temperature ramp
%experiment. Visualization of detrended traces and first derivative is
%intended to pick out activity of interneurons such as AIY

clear all; close all
disp ('Choose image acquisition file')
[image folder] = uigetfile('*.*');
imagefile = strcat(folder, image);

disp ('Choose image to select ROIs')
[file2, folder2]= uigetfile('*.*');
StanDev = strcat(folder2, file2);

nframes = input('How many frames is the image file? ');
freq = input('What is the frequency of acquisition in Hz? ');
maxintensity = input('Choose maximum pixel intensity for display');
ROI_choice_frame = input('Choose frame to pick ROIs');


M=multitiff2M_16bit(imagefile,1:nframes);
% Mr = registerfilm_ROI(M(:,:,1:nframes),1);
% flattenM = mean(M,3); 
x_dim = length(M(1,:,1));
y_dim = length(M(:,1,1));

figure;I=imread(StanDev,ROI_choice_frame);imshow(I,[0,maxintensity]);set(gcf, 'Position', get(0, 'Screensize'));
disp('Choose left set of ROIs')
% call function getROIcell which allows user to select ROIs
rois = getROIcell;
close all
    
% generate mask from rois, calculate dff
if isempty(rois);
    dff = []
    raw = []
    Mask = []
    int = [];
end

for i=1:size(rois,2); 
    Mask{i} = roipoly(imread(imagefile),rois{i}(:,1),rois{i}(:,2));
end;

%function calc_dff_NH calculates delta F/F for chosen ROIs
[dff, raw] = calc_dff_NH(M, Mask, size(rois,2));

%some of this is old code, but adj allows fit for photobleaching
%subtraction
% base = [];
% raw = [];
% int = [];
% adj = [];
% 
% nframes=size(dff,1);
% freq = 1;
% for i=1:size(dff,2); 
%     
%     [min_v, ind_v] = lmin(dff(:,i),1);
%     [fit_v, gof_v, out_v] = fit(ind_v', min_v', 'poly2', 'Normalize', 'on');
%     base (:,i) = feval(fit_v, [1:nframes])';
%     adj(:,i) = dff(:,i)-base(:,i);
%     int(:,i) = trapz(adj(:,i));
%     int(:,i) = int(:,i)/nframes*freq*60;
%     
% end

% Generate png/fig for mask of rois
figure; imshow(I,[0,maxintensity]);
for i=1:length(rois);
    patch(rois{1,i}(:,1),rois{1,i}(:,2),'m','FaceAlpha',0.3);
    text(rois{1,i}(1,1),rois{1,i}(1,2),num2str(i),'Color','g');
    hold on;
end

title('ROIs','FontSize', 18);
rois = strcat(folder, 'ROIs');
saveas(gcf, rois,'fig'); 
saveas(gcf, rois,'png');
close

%need this for Harry's script
mypath = fullfile(folder);
fname = image;

%Evidently this opens the txt file from labview, extracts temp vs time
%txt file needs to have the same base name as tif file
Lab_View_File = load(strcat(mypath,fname(1:end-3),'txt'));

Temperat = Lab_View_File(:,2);
Good_IND = [1:length(dff)];
Temperature = interp1(1:length(Temperat),Temperat,Good_IND);

Image_Time = length(dff);

resampletemp = Temperature;
pathname = mypath;

average_trace = mean(dff,2);

%If you want to average temperature traces uncomment below
%average_T = mean(Temperature);

figure('position', [200 0 2000 1200]); hold on;

for i=1:size(dff,2)
    
    yyaxis left; ylim([min(dff(:)), (max(dff(:)))]); F_trace = plot(dff(1:nframes,i),'-k'); F_trace.Color(4)=0.3;
    plot(average_trace(1:nframes),'-k','LineWidth',2);
    ylabel('\deltaF/F (%)');
    yyaxis right; ylim([(min(Temperature)-1), (max(Temperature)+1)]); 
    plot(Temperature(1:nframes),'-g','LineWidth',2);
    %or use this for average-
    %yyaxis right; ylim([(min(average_T)-1), (max(average_T)+1)]); 
    %plot(average_T,'-g','LineWidth',2);
    ylabel('Temperature');
end

vccd = strcat(folder, '1st_deriv');
saveas(gcf,vccd,'fig');
saveas(gcf,vccd,'png');

delta_F = strcat(folder, 'delta_F');
saveas(gcf, delta_F,'fig');

%Currently set up for 481 frame videos acquired at 2 Hz
%change "dff_cutoff" if you want to visualize a subsection of the trace,
%such as "Peri-Tc"
 
dff_cutoff = dff(40:481,:);
dff_cutoff_detrend = detrend(dff_cutoff,4);
 

fig6=figure; 

pos = get(fig6,'position');

set(fig6,'position',[pos(1:2)/4 pos(3:4)*2])

for i=1:size(dff_cutoff_detrend,2); 

    subplot((round(size(dff_cutoff_detrend,2)/2)),2,i);plot(dff_cutoff_detrend(:,i),'k'); hold on; 

    plot(diff(dff_cutoff_detrend(:,i)),'r');

    title(['ROI', num2str(i)]);%axis ([0 160 -10 30 ]);

    hold off; 

end;

vccd = strcat(folder, 'detrend');
saveas(gcf,vccd,'fig');
saveas(gcf,vccd,'png');

%delta_F = strcat(folder, 'delta_F');
%saveas(gcf, delta_F,'fig');

save analysis

