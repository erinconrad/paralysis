function method = load_paralysis_data(folder,file_name)


file_path = [folder,file_name];

[~, sheetNames] = xlsfinfo(file_path);
numSheets = length(sheetNames);


%% Get the raw data
cc_tbl = readtable(file_path, 'Sheet', 2);
cc_id = cc_tbl.record_id;
ep_true = cc_tbl.epileptiform_post_report;
art = cc_tbl.artifact_pre;

tbl = readtable(file_path,'Sheet',3);
pre_id = tbl.record_id;
pre = [tbl.R1_EC,tbl.R2_MG,tbl.R4_RR,tbl.R5_CE];
tbl = readtable(file_path,'Sheet',4);
pre_conf = [tbl.R1_EC,tbl.R2_MG,tbl.R4_RR,tbl.R5_CE];

tbl = readtable(file_path,'Sheet',5);
ar_id = tbl.record_id;
ar = [tbl.R1_EC,tbl.R2_MG,tbl.R4_RR,tbl.R5_CE];
tbl = readtable(file_path,'Sheet',6);
ar_conf = [tbl.R1_EC,tbl.R2_MG,tbl.R4_RR,tbl.R5_CE];

tbl = readtable(file_path,'Sheet',7);
par_id = tbl.record_id;
par = [tbl.R1_EC,tbl.R2_MG,tbl.R4_RR,tbl.R5_CE];
tbl = readtable(file_path,'Sheet',8);
par_conf = [tbl.R1_EC,tbl.R2_MG,tbl.R4_RR,tbl.R5_CE];

%% Sort by id to align everything
[cc_id,I] = sort(cc_id);
ep_true = ep_true(I);
art = art(I);

[pre_id,I] = sort(pre_id);
pre = pre(I,:);
pre_conf = pre_conf(I,:);

[ar_id,I] = sort(ar_id);
ar = ar(I,:);
ar_conf = ar_conf(I,:);

[par_id,I] = sort(par_id);
par = par(I,:);
par_conf = par_conf(I,:);

if ~isequal(cc_id,pre_id,ar_id,par_id)
    error('IDs do not align');
end

%% Put into organized structure
method(1).name = 'pre';
method(1).dc = logical(pre);
method(1).conf = logical(pre_conf);
method(1).ids = pre_id;

method(2).name = 'ar';
method(2).dc = logical(ar);
method(2).conf = logical(ar_conf);
method(2).ids = ar_id;

method(3).name = 'par';
method(3).dc = logical(par);
method(3).conf = logical(par_conf);
method(3).ids = par_id;

method(4).name = 'true';
method(4).dc = logical(ep_true);
method(4).ids = cc_id;
method(4).artifact = art;

end