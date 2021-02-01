function p_val = two_sample_bootstrap(all_reads,truth,artifact,all_acc,nb)

%{
I sample the EEGs with replacement to generate data for a confidence
interval
%}


twitch = artifact == 1;
emg = artifact == 0;

ntwitch = sum(twitch);
nemg = sum(emg);

diff_kappa_ar_boot = nan(nb,1);
diff_acc_ar_boot = nan(nb,1);

for ib = 1:nb
    
    twitch_acc = all_acc(:,twitch,:);
    emg_acc = all_acc(:,emg,:);
    
    %% Separately sample the twitch and emg eegs with replacement
    which_twitch = randi(ntwitch,ntwitch,1); % resample within the twitch eegs
    which_emg = randi(nemg,nemg,1);
    
    fake_twitch_acc = sum(sum(twitch_acc(:,which_twitch,:),2),1);
    fake_emg_acc = sum(sum(emg_acc(:,which_emg,:),2),1);
    
    %% take the difference between pre-AR and post-AR
    pre_post_twitch_acc = fake_twitch_acc(2) - fake_twitch_acc(1);
    pre_post_emg_acc = fake_emg_acc(2) - fake_emg_acc(1);
    
    %% Now get difference between emg and twitch (is improvement in twitch better than emg)?
    diff_pre_post_acc = pre_post_twitch_acc-pre_post_emg_acc;
    diff_acc_ar_boot(ib) = diff_pre_post_acc;
    %{
    twitch_reads = all_reads(:,twitch,:);
    emg_reads = all_reads(:,emg,:);
    
    
    twitch_truth = truth(twitch);
    emg_truth = truth(emg);
    
    %% Separately sample the twitch and emg eegs with replacement
    which_twitch = randi(ntwitch,ntwitch,1); % resample within the twitch eegs
    which_emg = randi(nemg,nemg,1);
    
    %% Get fake reads assuming these eegs
    fake_twitch = twitch_reads(:,which_twitch,:);
    fake_emg = emg_reads(:,which_emg,:);
    fake_truth_twitch = twitch_truth(which_twitch);
    fake_truth_emg = emg_truth(which_emg);
    
    
    
    %% Get kappas
    kappa_twitch = fleiss_kappa(fake_twitch);
    kappa_emg = fleiss_kappa(fake_emg);
    
    kappa_ar_diff = kappa_emg(2) - kappa_twitch(2); % is AR kappa better for emg than twitch?
    diff_kappa_ar_boot(ib) = kappa_ar_diff;
    %}
    
end

%% Sort the bootstrap accuracy differences
diff_acc_ar_boot = sort(diff_acc_ar_boot);
% Calculate how many are zero or less
num_zero_or_less_acc_ar = sum(diff_acc_ar_boot <= 0);
per_zero_or_less_acc_ar = num_zero_or_less_acc_ar/nb;

p_val.acc_ar = per_zero_or_less_acc_ar;

%{
%% Sort the bootstrap kappa differences
diff_kappa_ar_boot = sort(diff_kappa_ar_boot);
% Calculate how many are zero or less
num_zero_or_less_kappa_ar = sum(diff_kappa_ar_boot <= 0);
per_zero_or_less_kappa_ar = num_zero_or_less_kappa_ar/nb;

p_val.kappa_ar = per_zero_or_less_kappa_ar;
%}

end