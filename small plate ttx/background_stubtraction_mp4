clear all
[movieName, pathName] = uigetfile(['*.mp4']);
MovieFile = [pathName movieName];
[pathstr,name,ext] = fileparts(MovieFile); 
mov = VideoReader(MovieFile);
f = mov.NumberOfFrames;

disp('Reading images ...'); % create a background image by averaging all frames
pixValue = double(read(mov, 1)); % Read in the first frame of the video 
for i =  11:10:f % Read every 10 frames
    Frames = double(read(mov, i));
    pixValue = pixValue + Frames;
end

avgImage=uint8(pixValue/fix((f-11)/10+2));

disp('Subtracting background ...');% Subtract background from each image and make mp4 movie
bckSubV = VideoWriter([pathName filesep name '_bcksub'],'MPEG-4');
bckSubV.FrameRate = 1;
open(bckSubV);
for k = 1:f 
  this_frame = read(mov, k);
  bckSubImgs = bsxfun(@minus,this_frame,avgImage);
  writeVideo(bckSubV, bckSubImgs)
end
close(bckSubV);
disp('Done');

