function p_val = more_general_perm_test(all_reads,truth)

%{
all_reads is n_reviewers  x n_eegs x n_methods, 1 if epileptiform, 0 if not
truth is n_eegs
%}

nboot = 1e4;
n_reviewers = size(all_reads,1);
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

%% Get true accuracies
acc_baseline = sum(sum(accuracy(:,:,1)))/n_reads;
acc_ar = sum(sum(accuracy(:,:,2)))/n_reads;
diff_acc = acc_ar-acc_baseline;

%% Get true kappas
kappa_true = fleiss_kappa(all_reads);
kappa_true_diff = kappa_true(2) - kappa_true(1);

%% Get sensitivity and specificity
[sensitivity,specificity] = sens_and_spec(all_reads,truth);
sens_true_diff = sensitivity(2) - sensitivity(1);
spec_true_diff = specificity(2) - specificity(1);

%% Initialize permutation arrays
diff_acc_boot = nan(nboot,1);
diff_kappa_boot = nan(nboot,1);
diff_sens_boot = nan(nboot,1);
diff_spec_boot = nan(nboot,1);

baseline = all_reads(:,:,1);
ar = all_reads(:,:,2);

for ib = 1:nboot
    
    fake_baseline = baseline;
    fake_ar = ar;
    
    for j = 1:size(ar,2) % Loop over eegs
        r = randi([0 1]);
        if r == 1 % swap the whole column representing all reviewers for that read
            fake_baseline(:,j) = ar(:,j);
            fake_ar(:,j) = baseline(:,j);
        end
    end
    
    % Get accuracy
    fake_accuracy_baseline = is_read_correct(fake_baseline',truth);
    fake_accuracy_ar = is_read_correct(fake_ar',truth);
    acc_fake_baseline = sum(sum(fake_accuracy_baseline))/n_reads;
    acc_fake_ar = sum(sum(fake_accuracy_ar))/n_reads;
    diff_acc_boot(ib) = acc_fake_ar-acc_fake_baseline;
    
    % Get new kappa
    fake_reads(:,:,1) = fake_baseline;
    fake_reads(:,:,2) = fake_ar;
    fake_reads(:,:,3) = nan;
    kappa = fleiss_kappa(fake_reads);
    kappa_diff = kappa(2) - kappa(1);
    diff_kappa_boot(ib) = kappa_diff;
    
    % Get new sensitivity and specificity
    [fake_sens,fake_spec] = sens_and_spec(fake_reads,truth);
    fake_sens_diff = fake_sens(2)-fake_sens(1);
    fake_spec_diff = fake_spec(2)-fake_spec(1);
    diff_sens_boot(ib) = fake_sens_diff;
    diff_spec_boot(ib) = fake_spec_diff;
    
end

% Sort the bootstrap accuracy differences
diff_acc_boot = sort(diff_acc_boot);
n_more_accurate = sum(diff_acc_boot >= diff_acc);
p_val_acc = n_more_accurate/nboot;
if p_val_acc == 0 
    p_val_acc = 1/(nboot+1);
end

% Sort the bootstrap kappa differences
diff_kappa_boot = sort(diff_kappa_boot);
n_more_accurate = sum(diff_kappa_boot >= kappa_true_diff);
p_val_k = n_more_accurate/nboot;
if p_val_k == 0 
    p_val_k = 1/(nboot+1);
end

% Sort the bootstrap sens differences
diff_sens_boot = sort(diff_sens_boot);
n_more_accurate = sum(diff_sens_boot >= sens_true_diff);
p_val_sens = n_more_accurate/nboot;
if p_val_sens == 0 
    p_val_sens = 1/(nboot+1);
end

% Sort the bootstrap spec differences
diff_spec_boot = sort(diff_spec_boot);
n_more_accurate = sum(diff_spec_boot >= spec_true_diff);
p_val_spec = n_more_accurate/nboot;
if p_val_spec == 0 
    p_val_spec = 1/(nboot+1);
end

p_val.acc = p_val_acc;
p_val.kappa = p_val_k;
p_val.sens = p_val_sens;
p_val.spec = p_val_spec;

end



