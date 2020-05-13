%This scripts make the centers files that has 2100 or more images. If
%original data set has less than 2100 images, this will add blank images at
%the end.

clear all
[filename, fullpath] = uigetfile(['*.csv']);
Trackfile = [fullpath filesep filename];
[pathstr,name,ext] = fileparts(filename); 
Position = csvread(Trackfile, 5, 2);
Position(Position == 0) = -9;

T = [];
[H L] = size(Position);

%Add blank data to make total image number 2100
if (H<2100);
    Blank_H = 2100-H;
    Blank = zeros(Blank_H, L);
    Blank(:) = -9;
    Position = vertcat(Position, Blank);
    fprintf(2, '%d (less than 2100) images in %s and added blank images at the end\n', H, name);
    H = 2100;
end
    
for i = 1:H*2;
    
    if rem(i,2) == 1;
        R = Position((i+1)/2, 1:2:L);
       
    else
        R = Position(i/2, 2:2:L);
        
    end
    
    T = vertcat(T, R);
    i = i+1;
end
T_rev = T;
T_rev(2:2:length(T), :) = T(end:-2:1, :);
T_rev(1:2:length(T), :) = T(end-1:-2:1, :);
save([fullpath name '_centers.txt'], 'T_rev', '-ascii');

