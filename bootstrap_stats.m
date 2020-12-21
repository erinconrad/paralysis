function p_val = bootstrap_stats(all_reads,truth,all_conf,nb)

n_eegs = size(all_reads,2);
n_methods = size(all_reads,3);
n_reads = size(all_reads,1)*size(all_reads,2);

%% Get an accuracy matrix
% 1 if the read is correct, 0 if incorrect
accuracy = nan(size(all_reads));
for m = 1:n_methods
    reads = squeeze(all_reads(:,:,m));
    correct = is_read_correct(reads',truth);
    accuracy(:,:,m) = correct';
end


%% Initialize permutation arrays
diff_acc_boot = nan(nb,1);
diff_kappa_boot = nan(nb,1);
diff_sens_boot = nan(nb,1);
diff_spec_boot = nan(nb,1);
diff_conf_boot = nan(nb,1);

baseline = all_reads(:,:,1);
ar = all_reads(:,:,2);

baseline_conf = all_conf(:,:,1);
ar_conf = all_conf(:,:,2);

for ib = 1:nb
    
    %% Take n_eegs with replacement
    which_eegs = randi(n_eegs,n_eegs,1); % n_eegs x 1 matrix of which eegs to take
    
    %% Get fake reads assuming these eegs
    fake_baseline = baseline(:,which_eegs);
    fake_ar = ar(:,which_eegs);
    
    %% Get fake conf assuming these eegs
    fake_baseline_conf = baseline_conf(:,which_eegs);
    fake_ar_conf = ar_conf(:,which_eegs);
    
    %% Get new accuracies
    % Note that I need to update the truth array to contain whether the
    % corresponding new eegs have discharges or not
    fake_accuracy_baseline = is_read_correct(fake_baseline',truth(which_eegs));
    fake_accuracy_ar = is_read_correct(fake_ar',truth(which_eegs));
    
    
    acc_fake_baseline = sum(sum(fake_accuracy_baseline))/n_reads;
    acc_fake_ar = sum(sum(fake_accuracy_ar))/n_reads;
    diff_acc_boot(ib) = acc_fake_ar-acc_fake_baseline;
    
    %% Get new kappa
    fake_reads(:,:,1) = fake_baseline;
    fake_reads(:,:,2) = fake_ar;
    fake_reads(:,:,3) = nan;
    kappa = fleiss_kappa(fake_reads);
    kappa_diff = kappa(2) - kappa(1);
    diff_kappa_boot(ib) = kappa_diff;
    
    %% Get new sensitivity and specificity
    % Note that I need to update the truth array to contain whether the
    % corresponding new eegs have discharges or not
    [fake_sens,fake_spec] = sens_and_spec(fake_reads,truth(which_eegs));
    fake_sens_diff = fake_sens(2)-fake_sens(1);
    fake_spec_diff = fake_spec(2)-fake_spec(1);
    diff_sens_boot(ib) = fake_sens_diff;
    diff_spec_boot(ib) = fake_spec_diff;
    
    %% Get new conf
    diff_conf_boot(ib) = sum(sum(fake_baseline_conf))/n_reads - ...
        sum(sum(fake_ar_conf))/n_reads;
    
    
end

%% Sort the bootstrap accuracy differences
diff_acc_boot = sort(diff_acc_boot);
% Calculate how many are zero or less
num_zero_or_less_acc = sum(diff_acc_boot <= 0);
per_zero_or_less_acc = num_zero_or_less_acc/nb;

%% Sort the bootstrap kappa differences
diff_kappa_boot = sort(diff_kappa_boot);
% Calculate how many are zero or less
num_zero_or_less_kappa = sum(diff_kappa_boot <= 0);
per_zero_or_less_kappa = num_zero_or_less_kappa/nb;

%% Sort the bootstrap sens differences
diff_sens_boot = sort(diff_sens_boot);
% Calculate how many are zero or less
num_zero_or_less_sens = sum(diff_sens_boot <= 0);
per_zero_or_less_sens = num_zero_or_less_sens/nb;

%% Sort the bootstrap sens differences
diff_spec_boot = sort(diff_spec_boot);
% Calculate how many are zero or less
num_zero_or_less_spec = sum(diff_spec_boot <= 0);
per_zero_or_less_spec = num_zero_or_less_spec/nb;

%% Sort the bootstrap conf differences
diff_conf_boot = sort(diff_conf_boot);
% Calculate how many are zero or less
num_zero_or_less_conf = sum(diff_conf_boot <= 0);
per_zero_or_less_conf = num_zero_or_less_conf/nb;

p_val.acc = per_zero_or_less_acc;
p_val.kappa = per_zero_or_less_kappa;
p_val.sens = per_zero_or_less_sens;
p_val.spec = per_zero_or_less_spec;
p_val.conf = per_zero_or_less_conf;


end