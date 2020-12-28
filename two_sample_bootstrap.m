function p_val = two_sample_bootstrap(all_reads,truth,artifact,nb)

%{
I sample the EEGs with replacement to generate data for a confidence
interval
%}


twitch = artifact == 1;
emg = artifact == 0;

ntwitch = sum(twitch);
nemg = sum(emg);

diff_kappa_ar_boot = nan(nb,1);

for ib = 1:nb
    
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
    
end

%% Sort the bootstrap kappa differences
diff_kappa_ar_boot = sort(diff_kappa_ar_boot);
% Calculate how many are zero or less
num_zero_or_less_kappa_ar = sum(diff_kappa_ar_boot <= 0);
per_zero_or_less_kappa_ar = num_zero_or_less_kappa_ar/nb;

p_val.kappa_ar = per_zero_or_less_kappa_ar;

end