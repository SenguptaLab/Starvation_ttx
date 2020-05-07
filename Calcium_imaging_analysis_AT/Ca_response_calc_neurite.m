clear all; close all;

% load the image file
[fname,mypath]=uigetfile({'*.tif;*.mat;*.fig'},'Select the image file (tif) ...');
stack=tiffread2([mypath fname]);
[pathstr,filename,ext] = fileparts(fname); 

%extrack images from data structure
Total_Data = struct2cell(stack);
Images = Total_Data(7,1,:);
Images = cell2mat(Images);
L=size(Images);

%show projection of the stack and ROI
Image1(:,:) = max(Images,[],3);
colormap('default')
imagesc(Image1);

%Select background ROI
title('Select section of background...double click');
h = imrect;
BG_position = wait(h);
k = 1;
Analysis_folder = mkdir(mypath,'Analysis_files');

%Select ROI to analyze
while k>0  
      button = questdlg('Do you want to select ROI for analysis?','ROI selection','Yes','No','No');
      switch button
            case 'Yes'
                colormap('default')
                imagesc(Image1);
                title('Select the ROI for analysis...double click');
                h = imrect;
                ROI = wait(h);
                figure = imshow(Image1);
                hold on
                ROI_name = [filename,'_', num2str(k)];
                rectangle('Position', ROI,'EdgeColor','y');
                saveas(figure, [mypath 'Analysis_files/' ROI_name], 'jpg');
                close
                
                R_Neuron_X1 = uint16(ROI(1));
                R_Neuron_X2 = uint16(ROI(1) + ROI(3));
                R_Neuron_Y1 = uint16(ROI(2));
                R_Neuron_Y2 = uint16(ROI(2) + ROI(4));
                R_Neuron_BG_X1 = uint16(BG_position(1));
                R_Neuron_BG_X2 = uint16(BG_position(1) + BG_position(3));
                R_Neuron_BG_Y1 = uint16(BG_position(2));
                R_Neuron_BG_Y2 = uint16(BG_position(2) + BG_position(4));
                
                Neuron_Pixel_Size=20;
                R_Neuron_Im = Images(R_Neuron_Y1:R_Neuron_Y2, R_Neuron_X1:R_Neuron_X2, :);
                R_Neuron_Im = reshape(R_Neuron_Im,1,[],L(3));
                R_Neuron_Im = squeeze(R_Neuron_Im);
                R_Neuron_Im = R_Neuron_Im';
                R_Neuron_BG_Im = Images(R_Neuron_BG_Y1:R_Neuron_BG_Y2, R_Neuron_BG_X1:R_Neuron_BG_X2, :);
                R_Neuron_BG_Im = reshape(R_Neuron_BG_Im,1,[],L(3));
                R_Neuron_BG_Im = squeeze(R_Neuron_BG_Im);
                R_Neuron_BG = mean(R_Neuron_BG_Im);
                [R_Neuron, IX]= sort(R_Neuron_Im, 2,'descend');
                R_Neuron= R_Neuron(:,1:Neuron_Pixel_Size);
                M=size(R_Neuron_Im);
                R_Min_Col=R_Neuron(:,end);
                R_Neuron_loc=zeros(M(1),M(2));

                Right_Neuron_BG=zeros(M(1),M(2));

for i = 1:M(1)
   R_Neuron_loc(i,:)= R_Neuron_Im(i,:) >= R_Min_Col(i);
   Right_Neuron_BG(i,:)= R_Neuron_loc(i,:)*R_Neuron_BG(i);
end;

R_Neuron_Pixels=sum(R_Neuron_loc,2)';
Right_Neuron_BG=uint16(Right_Neuron_BG);
R_Neuron_loc=uint16(R_Neuron_loc);

R_Neuron_Im = R_Neuron_Im.*R_Neuron_loc;
Signal = R_Neuron_Im - Right_Neuron_BG;
Signal = sum(Signal,2)';
Norm_Signal=Signal./max(Signal);
Image_Time= find(Norm_Signal>=0);

Signal=Signal(Image_Time);

Right_Neuron_BG = sum(Right_Neuron_BG,2)'./R_Neuron_Pixels;
Right_Neuron_BG = Right_Neuron_BG(Image_Time);
 
Lab_View_File = load(strcat(mypath,fname(1:end-3),'txt'));
Temperat = Lab_View_File(:,2);
Temperature = interp1(1:length(Temperat),Temperat,Image_Time);
%Image_Time = Good_IND;

baseline=mean(Signal(1:10));
delta_F=((Signal-baseline)/baseline)*100;

[AX,H1,H2]=plotyy(Image_Time,delta_F, Image_Time,Temperature,'plot');
hold on
%[AX,H1,H2]=plotyy(Image_Time,FRET1, Image_Time,Temperature,'plot');
set(H1,'LineStyle','-','LineWidth',2,'color', 'b');
set(H2,'LineStyle','--','LineWidth',2,'color', 'g');
%ylim(axes_Delta_R,[-50 150])
ylabel('\deltaF/F (%)');
xlabel('Time (sec) ');
ylabel(AX(2),'Temperature \circC')

saveas(gcf,[mypath 'Analysis_files/' 'Delta_F_Figure_' ROI_name '.fig'],'fig');
saveas(gcf,[mypath 'Analysis_files/' 'Delta_F_Figure_' ROI_name '.jpg'],'jpg');
Mat_file_name = ['GCAMP_Analysis_' ROI_name '.mat'];
save ([mypath 'Analysis_files/' Mat_file_name],'baseline','BG_position','Signal','delta_F', 'Mat_file_name','h','i','Image_Time','mypath','ROI','Norm_Signal','R_Min_Col','R_Neuron_BG_X1','R_Neuron_BG_X2','R_Neuron_BG_Y1','R_Neuron_BG_Y2','R_Neuron_X1','R_Neuron_X2', 'R_Neuron_Y1','R_Neuron_Y2','Right_Neuron_BG','Temperat','Temperature','ROI_name')
        %clear AX BG_position button ext figure h H1 H2 Image1 Images IX
        %clear L_Neuron_loc R_Neuron_Im R_Neuron_loc IX L_IX L_Neuron_BG
        %clear L_Neuron_Pixels R_Neuron_BG R_Neuron_Pixels T_fileR_Min_Col M 
        %clear AX Delta_R_Figure H1 H2 M Neuron_Pixel_Size axes_Delta_R
         clf
%save ([mypath 'Analysis_files/' Mat_file_name], 'baseline',)
                k = k+1;
          case 'No'
             k = -1;
      end
end

close gcf