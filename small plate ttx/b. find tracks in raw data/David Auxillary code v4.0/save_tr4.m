function path_and_fn = save_tr4(data_path,pref,tr)
% save tr in data_path, 
% in a file name that starts with the string in pref and has a time-stamp

if data_path(end)=='\' % get rid of trailing '\'
    dp = data_path(1:(end-1));
else 
    dp = data_path; 
end; 

tmp = strfind(dp,'\'); % find last folder in data_path (dp) for file_name
if isempty(tmp)
    fn = [pref dp '_'];
else
	fn = [pref dp((tmp(end)+1):end) '_'];
end;
tmp = datestr(now,31); % add time-stamp to file_name
tmp = ['D' tmp(1:end-3)]; % get rid of seconds
tmp = regexprep(tmp,{':'},'_');
tmp = regexprep(tmp,{' '},'_H');
fn = [fn tmp '.mat']; % add extension to file-name

path_and_fn = [dp '\' fn];
save(path_and_fn, 'tr'); % save tr in file_name in data_path

return;
