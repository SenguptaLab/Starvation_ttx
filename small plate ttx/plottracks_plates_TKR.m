clear; close all;
[fle,pth] = uigetfile('/Volumes/home/bev/behavior/*.mat');
load([pth fle]);
 [a b] = fileparts(fileparts(pth));
pltid = [data_summary.tr.plate_id];
flid = {data_summary.tr.file_name_id};
color_by_time = 1;


if(color_by_time)
clrs = hsv(3813);
end


for p = 1:max([data_summary.tr.plate_id])
    figure('Color',[1 1 1]);
    hold on
    for i = 1:length(data_summary.tr)
        if(data_summary.tr(i).plate_id == p)
            
        if(~color_by_time)    
        clrs = hsv(length(data_summary.tr(i).x));
        end
        temps = interp1(linspace(data_summary.tr(i).min_edge,data_summary.tr(i).max_edge),linspace(data_summary.tr(i).min_T,data_summary.tr(i).max_T),data_summary.tr(i).x);
        %hold on
        for j = 2:length(temps)
            if(~color_by_time)    
            line([temps(j-1) temps(j)],[data_summary.tr(i).y(j-1) data_summary.tr(i).y(j)],'Color',clrs(j,:))
            else
            line([temps(j-1) temps(j)],[data_summary.tr(i).y(j-1) data_summary.tr(i).y(j)],'Color',clrs(data_summary.tr(i).f(j),:))
            end
            
            % pause(.1)
            if(j == length(temps))
               scatter(temps(j),data_summary.tr(i).y(j),30,[0 0 0],'*')
            end
            if(j == 2)
               scatter(temps(j-1),data_summary.tr(i).y(j-1),30,clrs(j,:),'o')
            end
        end
        %hold off
        end
      
    end
    gray = 0.7;
    set(gca,'Color',[gray gray gray],'YTick',[],'YColor',[gray gray gray]);
    hold off
      f = flid(pltid == p);
        cr =  data_summary.plates.I_cr(p);
      title([b f(1)  num2str(cr)],'Interpreter','None');
        if(color_by_time)
            colormap(clrs);
            set(gca,'CLim',[0 size(clrs,1)]);
            set(gca,'XLim',[18 23]);
            cb = colorbar('YLim',[0 4200],'YTick',linspace(0,4200,8),'YTickLabel', [' 0';'10';'20';'30';'40';'50';'60';'70']);
        end
    
end


