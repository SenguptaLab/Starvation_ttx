function [] = make_plot_contribution_of_runs_to_Icr(ds)

[N1,centers1,N3,centers3] = bin_runs_up_down4(ds);

[dur,ang,indx_dn,indx_up] = group_runs_up_down4(ds);

figure; 
clf;
bin_limits = [0 (centers1(2:end)+centers1(1:(end-1)))/2 999]; 
for k = 1:(length(bin_limits)-1)
    runs_in_bin_k = find(dur>=bin_limits(k) & dur<bin_limits(k+1));
    dn_bin_k = intersect(runs_in_bin_k,indx_dn);
    up_bin_k = intersect(runs_in_bin_k,indx_up);
    delta_by_bins(k) = sum(dur(dn_bin_k))-sum(dur(up_bin_k));
    runs_up_present(k) = sum(dur(up_bin_k))>0;
end;
plot(centers1,delta_by_bins/sum(delta_by_bins),'x:','markersize',8);
hold on;
plot(centers1(runs_up_present),delta_by_bins(runs_up_present)/sum(delta_by_bins),'or','markersize',8);
hold off;
title(['fraction contribution of each bin to I_{cr}']);
xlabel(['Mean duration of binned runs (sec)']);
ylabel(['Fraction of contribution to cryophilic bias']);
return;
