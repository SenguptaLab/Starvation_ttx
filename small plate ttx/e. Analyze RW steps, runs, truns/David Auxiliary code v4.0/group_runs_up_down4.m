function [dur,ang,indx_dn,indx_up] = group_runs_up_down4(ds)

global ANG_FOR_RW

dur = [];
ang = [];
for k = 1:length(ds.runs.time)    
    dur = [dur; ds.runs.time{k}];
    ang = [ang; ds.runs.angle{k}];
end;

% group symmetric angles so all are in range [0,pi]
ang(ang>pi) = 2*pi-ang(ang>pi); % there shouldn't be any angels >pi but..   
ang = abs(ang); % in terms of up/down gradient alpha is equivalent to -alpha 

indx_dn = find(ang>(pi-ANG_FOR_RW) & ang<=pi);
%indx_vr = find(ang<(pi/2+ANG_FOR_RW) & ang>pi/2-ANG_FOR_RW);
indx_up = find(ang<ANG_FOR_RW & ang>=0);

return;
