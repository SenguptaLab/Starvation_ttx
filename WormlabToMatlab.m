clear all
[filename, fullpath] = uigetfile(['*.csv']);
Trackfile = [fullpath filesep filename];
[pathstr,name,ext] = fileparts(filename); 
Position = csvread(Trackfile, 5, 2);
Position(Position == 0) = -9;

T = [];
[H L] = size(Position);
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

