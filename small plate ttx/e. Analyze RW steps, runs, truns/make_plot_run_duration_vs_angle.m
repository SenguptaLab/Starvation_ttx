function [] = make_plot_run_duration_vs_angle(ds)

global ANG_BIN

SHAPE='o';  % shape in plot

[dur,ang,indx_dn,indx_up] = group_runs_up_down4(ds);

edges = [0:ANG_BIN:pi];
cents = edges(1) + (ANG_BIN/2:ANG_BIN:edges(end));
[n,inds] = histc(ang,edges);

u = 1:length(edges);
for i = 1:length(u)-1;
    choose = find(inds==u(i));
    mdur(i) = mean(dur(choose));
    sdur(i) = std(dur(choose))/sqrt(length(choose));
end;

figure; 
clf;
hold on;
h = plot(cents*180/pi,mdur,[SHAPE 'k-']); 
for i = 1:length(cents)
    h = line(180/pi*[cents(i),cents(i)],[mdur(i)-sdur(i),mdur(i)+sdur(i)]); 
    set(h, 'color', 'k');
end;
ylabel('Mean run duration (sec)');
xlabel(['Angle of run (deg)']);

%%% fits to various models (Alex/Damon code)
%[pf,sf]  = lsq('modelcos', [1 1], cents, mdur, sdur, 0);     % cos(.) fit
%[pf2,sf2]= lsq('modelsgncos', [1 1], cents, mdur, sdur, 0);  % sign(cos(.)) fit
%[pl,sl]  = lsq('modellinear', [1 1], cents, mdur, sdur, 0);  % line fit
%[pc,sc]  = lsq('modelconst', [30], cents, mdur, sdur, 0);    % const fit
%title(['P_{const}=' num2str(sc.p_fit) ', P_{cos}= ' num2str(sf.p_fit) ...
%        ', P_{sgn(cos)}= ' num2str(sf2.p_fit) ', P_{lin}= ' num2str(sl.p_fit)]);
%set(gca,'ylim',[10,40]);

return;