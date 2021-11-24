function p_val = two_sample_bootstrap(all_reads,truth,artifact,all_acc,nb)

%{
I sample the EEGs with replacement to generate data for a confidence
interval
%}



emg = artifact == 1;
emg_no = artifact == 0;

nemg = sum(emg);
nemg_no = sum(emg_no);

%% Separate EEGs with lots and little emg
emg_acc = all_acc(:,emg,:);
emg_no_acc = all_acc(:,emg_no,:);
emg_reads = all_reads(:,emg,:);
emg_no_reads = all_reads(:,emg_no,:);
emg_truth = truth(emg);
emg_no_truth = truth(emg_no);

diff_kappa_ar_boot = nan(nb,1);
diff_acc_ar_boot = nan(nb,1);

for ib = 1:nb
     
    %% Separately sample the lots and little emg eegs with replacement
    which_emg = randi(nemg,nemg,1); % resample within the twitch eegs
    which_emg_no = randi(nemg_no,nemg_no,1);
    
    %% Get fake reads assuming these eegs
    fake_baseline_emg = emg_reads(:,which_emg,1);
    fake_ar_emg = emg_reads(:,which_emg,2);
    fake_par_emg = emg_reads(:,which_emg,3);
    fake_truth_emg = emg_truth(which_emg); % also get the truth about dc vs no assuming these eegs
    
    fake_baseline_no_emg = emg_no_reads(:,which_emg_no,1);
    fake_ar_no_emg = emg_no_reads(:,which_emg_no,2);
    fake_par_no_emg = emg_no_reads(:,which_emg_no,3);
    fake_truth_no_emg = emg_no_truth(which_emg_no); % also get the truth about dc vs no assuming these eegs
    
    %% Make new fake lots and little emg EEGs from the resampled EEGs
    fake_emg_acc = sum(sum(emg_acc(:,which_emg,:),2),1);
    fake_emg_no_acc = sum(sum(emg_no_acc(:,which_emg_no,:),2),1);
    
    %% take the difference between pre-AR and post-AR 
    pre_post_emg_acc = fake_emg_acc(2) - fake_emg_acc(1);
    pre_post_emg_no_acc = fake_emg_no_acc(2) - fake_emg_no_acc(1);
    
    
    %% Now get difference between lots and little emg (is improvement in lots emg better than little emg)?
    diff_pre_post_acc = pre_post_emg_acc-pre_post_emg_no_acc;
    diff_acc_ar_boot(ib) = diff_pre_post_acc;
    
    %% Get kappas
    fake_reads_emg(:,:,1) = fake_baseline_emg;
    fake_reads_emg(:,:,2) = fake_ar_emg;
    fake_reads_emg(:,:,3) = fake_par_emg;
    kappa_emg = fleiss_kappa(fake_reads_emg);
    pre_post_emg_kappa = kappa_emg(2)-kappa_emg(1);
    
    fake_reads_no_emg(:,:,1) = fake_baseline_no_emg;
    fake_reads_no_emg(:,:,2) = fake_ar_no_emg;
    fake_reads_no_emg(:,:,3) = fake_par_no_emg;
    kappa_no_emg = fleiss_kappa(fake_reads_no_emg);
    pre_post_no_emg_kappa = kappa_no_emg(2)-kappa_no_emg(1);
    
    diff_kappa_ar_boot(ib) = [pre_post_emg_kappa-pre_post_no_emg_kappa];
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

p_val.acc_ar = 2*per_zero_or_less_acc_ar;

%
%% Sort the bootstrap kappa differences
diff_kappa_ar_boot = sort(diff_kappa_ar_boot);
% Calculate how many are zero or less
num_zero_or_less_kappa_ar = sum(diff_kappa_ar_boot <= 0);
per_zero_or_less_kappa_ar = num_zero_or_less_kappa_ar/nb;

p_val.kappa_ar = per_zero_or_less_kappa_ar;

if 0
    plot(diff_kappa_ar_boot,'o')
    hold on
    plot(xlim,[0 0],'k--')
    title(sprintf('%1.3f',pval.kappa_ar))
end
%}

end