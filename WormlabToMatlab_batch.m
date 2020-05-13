clear all;
close all;
dir_of_files = uigetdir();
addpath(dir_of_files);
D = dir([dir_of_files filesep '*.csv']);
for i = 1:length(D);
    Filename = strsplit(D(i).name, '_');
    tf = strcmp(Filename(end), 'cal.csv');
    if tf == 1;
       i = i + 1;
    
    else
    fprintf('converting %s\n', D(i).name);   
    Position = csvread(D(i).name, 5, 2);
    Position(Position == 0) = -9;
     
    T = [];
    [H L] = size(Position);
    
    if (H<2100);
        
     
    
    for k = 1:H*2;
        if rem(k,2) == 1;
        R = Position((k+1)/2, 1:2:L);
        else
        R = Position(k/2, 2:2:L);
        end
        
        
        T = vertcat(T, R);
        k = k+1;
    
    end
    
    T_rev = T;
    T_rev(2:2:length(T), :) = T(end:-2:1, :);
    T_rev(1:2:length(T), :) = T(end-1:-2:1, :);
    save([dir_of_files filesep D(i).name(1:end-4) '_centers.txt'], 'T_rev', '-ascii');

    i = i + 1;
    end

        
    end

