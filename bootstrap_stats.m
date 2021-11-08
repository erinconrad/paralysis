function p_val = bootstrap_stats(all_reads,truth,all_conf,nb,two_stage)

%{
This is the main statistical test for the paralysis project.

I sample the EEGs with replacement to generate data for a confidence
interval.
%}

%% Get numbers of eegs, raters, and total reads
n_eegs = size(all_reads,2);
n_raters = size(all_reads,1);
n_reads = size(all_reads,1)*size(all_reads,2);

%% Initialize permutation arrays
diff_acc_boot = nan(nb,2); % first is AR-baseline, second is paralysis-baseline
diff_kappa_boot = nan(nb,2);
diff_sens_boot = nan(nb,2);
diff_spec_boot = nan(nb,2);
diff_conf_boot = nan(nb,2);

% These are the reads of dc vs no dc
baseline = all_reads(:,:,1);
ar = all_reads(:,:,2);
par = all_reads(:,:,3);

% The confidence
baseline_conf = all_conf(:,:,1);
ar_conf = all_conf(:,:,2);
par_conf = all_conf(:,:,3);

% Loop over nb bootstrap samples
for ib = 1:nb
    
    %% Take n_eegs with replacement
    which_eegs = randi(n_eegs,n_eegs,1); % n_eegs x 1 matrix of which eegs to take
    
    %% Get fake reads assuming these eegs
    fake_baseline = baseline(:,which_eegs);
    fake_ar = ar(:,which_eegs);
    fake_par = par(:,which_eegs);
    fake_truth = truth(which_eegs); % also get the truth about dc vs no assuming these eegs
    
    % I am not doing this
    if two_stage == 1
        %% Now, randomly sample with replacement n_reads
        % This is the 2nd stage of two-stage bootstrapping
        which_raters = randi(n_raters,n_raters,1);
        fake_baseline = fake_baseline(which_raters,:);
        fake_ar = fake_ar(which_raters,:);
        fake_par = fake_par(which_raters,:);
    end
    
    %% Get fake conf assuming these eegs
    fake_baseline_conf = baseline_conf(:,which_eegs);
    fake_ar_conf = ar_conf(:,which_eegs);
    fake_par_conf = par_conf(:,which_eegs);
    
    %% Get new accuracies
    % Note that I need to update the truth array to contain whether the
    % corresponding new eegs have discharges or not
    fake_accuracy_baseline = is_read_correct(fake_baseline',fake_truth);
    fake_accuracy_ar = is_read_correct(fake_ar',fake_truth);
    fake_accuracy_par = is_read_correct(fake_par',fake_truth);
    
    % Get the total proportion of correct eegs
    acc_fake_baseline = sum(sum(fake_accuracy_baseline))/n_reads;
    acc_fake_ar = sum(sum(fake_accuracy_ar))/n_reads;
    acc_fake_par = sum(sum(fake_accuracy_par))/n_reads;
    
    % Take the difference between the AR and the baseline, and between the
    % paralyzed and the baseline
    diff_acc_boot(ib,:) = [acc_fake_ar-acc_fake_baseline,...
        acc_fake_par - acc_fake_baseline];
    
    %% Get new kappa
    fake_reads(:,:,1) = fake_baseline;
    fake_reads(:,:,2) = fake_ar;
    fake_reads(:,:,3) = fake_par;
    kappa = fleiss_kappa(fake_reads);
    diff_kappa_boot(ib,:) = [kappa(2)-kappa(1),kappa(3)-kappa(1)];
    
    %% Get new sensitivity and specificity
    % Note that I need to update the truth array to contain whether the
    % corresponding new eegs have discharges or not
    [fake_sens,fake_spec] = sens_and_spec(fake_reads,fake_truth);
    diff_sens_boot(ib,:) = [fake_sens(2)-fake_sens(1),fake_sens(3)-fake_sens(1)];
    diff_spec_boot(ib,:) = [fake_spec(2)-fake_spec(1),fake_spec(3)-fake_spec(1)];
    
    %% Get new conf
    diff_conf_boot(ib,:) = [sum(sum(fake_ar_conf))/n_reads - ...
        sum(sum(fake_baseline_conf))/n_reads,...
        sum(sum(fake_par_conf))/n_reads - ...
        sum(sum(fake_baseline_conf))/n_reads];
    
    
end

%% Sort the bootstrap accuracy differences
diff_acc_boot = sort(diff_acc_boot,1);

% Calculate how many are zero or less. Note that this is a one-sided
% p-value. 
num_zero_or_less_acc = sum(diff_acc_boot <= 0,1);
per_zero_or_less_acc = num_zero_or_less_acc/nb;

if 0
    figure
    plot(diff_acc_boot(:,1),'o')
    hold on
    plot(xlim,[0 0],'k--')
    title(sprintf('p = %1.3f',per_zero_or_less_acc(1)))
end




%% Sort the bootstrap kappa differences
diff_kappa_boot = sort(diff_kappa_boot,1);
% Calculate how many are zero or less
num_zero_or_less_kappa = sum(diff_kappa_boot <= 0,1);
per_zero_or_less_kappa = num_zero_or_less_kappa/nb;

%% Sort the bootstrap sens differences
diff_sens_boot = sort(diff_sens_boot,1);
% Calculate how many are zero or less
num_zero_or_less_sens = sum(diff_sens_boot <= 0,1);
per_zero_or_less_sens = num_zero_or_less_sens/nb;

%% Sort the bootstrap spec differences
diff_spec_boot = sort(diff_spec_boot,1);
% Calculate how many are zero or less
num_zero_or_less_spec = sum(diff_spec_boot <= 0,1);
per_zero_or_less_spec = num_zero_or_less_spec/nb;

%% Sort the bootstrap conf differences
diff_conf_boot = sort(diff_conf_boot,1);
% Calculate how many are zero or less
num_zero_or_less_conf = sum(diff_conf_boot <= 0,1);
per_zero_or_less_conf = num_zero_or_less_conf/nb;

p_val.acc = per_zero_or_less_acc;
p_val.kappa = per_zero_or_less_kappa;
p_val.sens = per_zero_or_less_sens;
p_val.spec = per_zero_or_less_spec;
p_val.conf = per_zero_or_less_conf;


end