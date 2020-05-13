function Plot_Tracks_and_turns(fgn, trck)

if length(trck)>9
    INDX = unique(round(1+rand(1,20)*(length(trck)-1)));
    if length(INDX)>9
        INDX = INDX(1:9);
    end;
else
    INDX = 1:length(trck);
end;

sp = 0;
for i= INDX % 1:9 
    sp = sp+1;
    figure(fgn); subplot(3,3,sp);
    %figure(fgn+sp);
    hold on;
    if i<=length(trck)
        tr = trck(i); 
        for k = 1:length(tr.x) % plot original track
            sc = k/length(tr.x);
            plot(tr.x(k), tr.y(k), '.', 'color', [0 0 sc]);
        end;
        for m = 1:tr.run_num % show "runs" as straight lines (purple)
        	x1 = tr.x(tr.run_indx(m,1));
            y1 = tr.y(tr.run_indx(m,1));
            x2 = tr.x(tr.run_indx(m,2));
            y2 = tr.y(tr.run_indx(m,2));
            lh = plot([x1,x2],[y1,y2],'-','color',[0.5 0 0.5],'linewidth',6);
        end;

        for m = 1:2:tr.seg_num % show odd segments as straight lines (green)
        	x1 = tr.x(tr.seg_indx(m,1));
            y1 = tr.y(tr.seg_indx(m,1));
            x2 = tr.x(tr.seg_indx(m,2));
            y2 = tr.y(tr.seg_indx(m,2));
            lh = plot([x1,x2],[y1,y2],'-','color',[0 0.8 0],'linewidth',2);
        end;
        for m = 2:2:tr.seg_num % show even segments as straight lines (yellow)
        	x1 = tr.x(tr.seg_indx(m,1));
            y1 = tr.y(tr.seg_indx(m,1));
            x2 = tr.x(tr.seg_indx(m,2));
            y2 = tr.y(tr.seg_indx(m,2));
            lh = plot([x1,x2],[y1,y2],'-','color',[1 1 0],'linewidth',2);
        end;
        % show approxmate turns locations as circles (red)
        plot(tr.x(tr.turn_indx),tr.y(tr.turn_indx),'o','color',[1 0 0]);
        title(([num2str(tr.run_num) ' runs']));
    end;
    hold off;
end;

return;