clear; close all;
[fle,pth] = uigetfile('/Volumes/home/bev/behavior/*.mat');
load([pth fle]);
 [a b] = fileparts(fileparts(pth));
pltid = [data_summary.tr.plate_id];
flid = {data_summary.tr.file_name_id};

prompt = {'Enter temperature on the lower side:','Enter temperature on the higher side:'};
dlg_title = 'Temperature setting';
num_lines = 1;
answer = inputdlg(prompt,dlg_title,num_lines);
Low_Temp = str2double(answer{1});
High_Temp = str2double(answer{2});

clrs = magma(2100);   %Change this to jet for heatmap-like color

button = questdlg('Do you want to indicate beginning and end of tracks?','Track appearance','Yes','No','Yes');
      switch button
            case 'Yes'
                Indication = 1;
            case 'No'
                Indication = 0;
      end
      
              
for p = 1:max([data_summary.tr.plate_id]);
    figure %('Color',[1 1 1]);
    hold on
    if (Low_Temp == High_Temp);
        for k = 1:length(data_summary.tr);
            if (data_summary.tr(k).plate_id == p)
                point_number = length(data_summary.tr(k).x);
                for m = 2:point_number;
                   line([data_summary.tr(k).x(m-1) data_summary.tr(k).x(m)],[data_summary.tr(k).y(m-1) data_summary.tr(k).y(m)],'Color',clrs(data_summary.tr(k).f(m),:));
                   if (m == point_number)&&(Indication == 1);
                      scatter(data_summary.tr(k).x(m),data_summary.tr(k).y(m),30,[0 0 0],'*');
                   end;
                   if(m == 2)&&(Indication == 1);
                       scatter(data_summary.tr(k).x(m-1),data_summary.tr(k).y(m-1),30,clrs(m,:),'o');
                   end
                end
            end
        end
    else
 
        for i = 1:length(data_summary.tr)
            if(data_summary.tr(i).plate_id == p)
                temps = interp1(linspace(data_summary.tr(i).min_edge,data_summary.tr(i).max_edge),linspace(data_summary.tr(i).min_T,data_summary.tr(i).max_T),data_summary.tr(i).x);
        
        for j = 2:length(temps)
            line([temps(j-1) temps(j)],[data_summary.tr(i).y(j-1) data_summary.tr(i).y(j)],'Color',clrs(data_summary.tr(i).f(j),:))
            
            if(j == length(temps))&&(Indication == 1)
               scatter(temps(j),data_summary.tr(i).y(j),30,[0 0 0],'*')
            end
            if(j == 2)&&(Indication == 1)
               scatter(temps(j-1),data_summary.tr(i).y(j-1),30,clrs(j,:),'o')
            end
        end
        
            end
        end
        
    end
    
    gray = 0.7;
    set(gca,'Color',[gray gray gray],'YTick',[],'YColor',[gray gray gray]);
    hold off
      f = flid(pltid == p);
        cr =  data_summary.plates.I_cr(p);
      title([b f(1)  num2str(cr)],'Interpreter','None');
      
            colormap(clrs);
            set(gca,'CLim',[0 size(clrs,1)]);
            if isequal(Low_Temp, High_Temp)==0 
            set(gca,'XLim',[Low_Temp High_Temp]);
            else
               set(gca,'XLim',[0 1280]);  %set as the width of picture in px
               set(gca,'YLim',[0 1024]);   %set as the height of picture in px
            end
            cb = colorbar('YLim',[0 2100],'YTick',linspace(0,2100,8),'YTickLabel', [' 0';' 5';'10';'15';'20';'25';'30';'35']);
        
end

 

