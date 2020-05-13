clear all;
close all;
c = clock;
path = ['C:\ThermoAssays\' num2str(c(2)) '.' num2str(c(1)) filesep];
fn = 1;
while ~isequal(fn,0) % 0 = output value of uigetfile if user clicks "cancel"
    [fn,path] = uigetfile([path '*.mat'],'Select mat file; or click "cancel"');
     if ~isequal(fn,0) % a file was chosen
         if  exist([path filesep 'BLC'], 'dir') ==0;
             mkdir(path,'BLC');
         else 
             disp('BLC file added to the folder already exists');
         end
     
        load([path fn]);
        set(0,'DefaultFigurePosition',[100 100 1000 1000]);
        fig  = plot(delta_F);
        disp(fn);
        title(fn,'Interpreter','none');
        ax =gca;
        axes(ax);
        [blx, bly] = getline(ax);
        [ublx, a, b] = unique(blx);
        ubly = bly(a);

        bsl = interp1(ublx,ubly,1:length(delta_F));
        bsl(isnan(bsl)) = 0;

        close

        BLC_raw_delta_F = (delta_F + abs(min(bly))) - (bsl +  abs(min(bly)));

        figure;
        hold on
        plot(BLC_raw_delta_F,'r');
        title('Corrected signal')
        pause(2)
        close

        BLC_file = (['blcrrt_', ROI_name]);
        save([path, 'BLC/', BLC_file]);

     end
end %if file is not chosen



         
        


