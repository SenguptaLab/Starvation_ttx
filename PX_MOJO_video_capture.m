%%Script acquires video at  1fps

vid = videoinput('pixelinkimaq', 1, 'MONO8_1280x1024');

src = getselectedsource(vid);
src.Exposure = 100;
src.Gain = '1.14';
src.GammaMode = 'on';
src.Gamma = 0.100001;
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
vid.LoggingMode = 'disk';
folder_name = uigetdir();
diskLogger = VideoWriter([folder_name filesep filename], 'MPEG-4');
v.Quality = 100;


diskLogger.FrameRate = 1; 
vid.DiskLogger = diskLogger;
preview(vid);

button = questdlg('Click OK to start video!','Video Start Box','OK','cancel','cancel');
switch button
    case 'OK'
        run = true;
        display('Recording started.')
        start(vid);
            case 'cancel'
        run = false;
        display('Recording canceled.')
end
%%