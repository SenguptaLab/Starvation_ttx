%%Script acquires video at  1fps

vid = videoinput('pixelinkimaq', 1, 'MONO8_1280x1024');

src = getselectedsource(vid);
src.Exposure = 100;
src.Gain = '2.48';
src.GammaMode = 'on';
src.Gamma = 1.0085;
src.ActualFrameRate = 1;
src.FrameRate = 1;
vid.ReturnedColorspace = 'grayscale';

prompt = {'Enter Videolength(frames):','Enter Filename:'};
dlg_title = 'Video Specifications';
num_lines = 1;
def = {'2100',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
filename = answer{2};

vid.FramesPerTrigger = str2double(answer{1}); %3600 = 30min at 2fps
% vid.FrameGrabInterval = 7;
vid.LoggingMode = 'disk&memory';
imaqmem(5000000000); %set maximum allowed memory usage
folder_name = uigetdir();
diskLogger = VideoWriter([folder_name filesep filename], 'MPEG-4');
%v.Quality = 100;


diskLogger.FrameRate = 1; 
vid.DiskLogger = diskLogger;
preview(vid);

button = questdlg('Click OK to start video!','Video Start Box','OK','cancel','cancel');
switch button
    case 'OK'
        run = true;
        display('Recording started.');
        start(vid)
        wait(vid, inf);
      
    case 'cancel'
        run = false;
        display('Recording canceled.')
end

%Background subtraction
disp('Reading images ...'); % create a background image by averaging all frames
f = getdata(vid);
Num = size(f, 4);
pixValue = double(f(:,:,:,1)); % Read in the first frame of the video 
for i = 11:10:Num % Read every 10 frames
    Frames = double(f(:,:,:,i));
    pixValue = pixValue + Frames;
end

avgImage=uint8(pixValue/fix((Num-11)/10+2));

disp('Subtracting background ...');% Subtract background from each image and make mp4 movie
bckSubV = VideoWriter([diskLogger.Path filesep filename '_bcksub'],'MPEG-4');
bckSubV.FrameRate = 1;
open(bckSubV);
for k = 1:Num 
  this_frame = f(:,:,:,k);
  bckSubImgs = bsxfun(@minus,this_frame,avgImage);
  writeVideo(bckSubV, bckSubImgs)
end
close(bckSubV);
disp('Done');


%%