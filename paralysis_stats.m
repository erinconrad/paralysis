clear

%{
I added some papers to paperpile suggesting that a single stage bootstrap
analysis is a reasonable way to handle clustered data.


%}

%% Parameters
two_stage = 0;
data_folder = '../data/';
file_name = 'Persyst Data Raw 02-04-21.xls';
results_folder = '../results/';

%% Get data
method = load_paralysis_data(data_folder,file_name);

%% Convert epileptiform/no-epileptiform to correct/incorrect
for i = 1:3 % Loop over methods
    method(i).correct = is_read_correct(method(i).dc,method(4).dc); % method(4) is report
end

%{
%% Concatenate into one long vector (treating all reads as independent)
for i = 1:3 % Loop over methods
    method(i).correct_vec = method(i).correct(:);
    method(i).conf_vec = method(i).conf(:);
end
%}

%% Convert baseline reads into one matrix
all_ep_or_no = nan(size(method(1).correct,2),size(method(1).correct,1),3);
% n_reviewers x n_eegs x n_methods
for i = 1:3 % loop over methods
    all_ep_or_no(:,:,i) = method(i).dc'; % transpose to convert to n_reviewers x n_eegs
end

%% Convert confidence into one matrix
all_conf = nan(size(method(1).correct,2),size(method(1).correct,1),3);
% n_reviewers x n_eegs x n_methods
for i = 1:3
    all_conf(:,:,i) = method(i).conf';
end

%% Convert accuracy into one matrix
% Convert into one single 3 dimensional matrix n_reviewers  x n_eegs x
% n_methods
all_acc = nan(size(method(1).correct,2),size(method(1).correct,1),3);
for i = 1:3 % loop over methods
    all_acc(:,:,i) = method(i).correct'; % transpose to convert to n_reviewers x n_eegs
end
n_eegs = size(all_ep_or_no,2);

%% Basic info
% How many patients are seizing?
fprintf('\n\n%d of %d patients had epileptiform discharges per clinical report.\n',...
    sum(method(4).dc),n_eegs);

%% Bootstrap statistics
pval = bootstrap_stats(all_ep_or_no,method(4).dc,all_conf,1e4,two_stage);

%% Compare percentages correct
n_reads = length(method(1).correct(:));
fprintf('\n\nAverage percent reads that are correct:\n');
fprintf('Baseline: %1.1f%%\n',sum(method(1).correct(:))/n_reads*100);
fprintf('AR: %1.1f%%\n',sum(method(2).correct(:))/n_reads*100);
fprintf('Paralysis: %1.1f%%\n',sum(method(3).correct(:))/n_reads*100);
fprintf('Accuracy p-value: %1.3f\n',pval.acc);


%% Compare kappa statistics
kappa = fleiss_kappa(all_ep_or_no);
fprintf('\n\nKappa statistics:\n');
fprintf('Baseline: %1.2f\n',kappa(1));
fprintf('AR: %1.2f\n',kappa(2));
fprintf('Paralysis: %1.2f\n',kappa(3));
fprintf('Kappa p-value: %1.3f\n',pval.kappa);

%% Compare sensitivity and specificity
[sensitivity,specificity] = sens_and_spec(all_ep_or_no,method(4).dc);
fprintf('\n\nSensitivity and specificity (respectively):\n');
fprintf('Baseline: %1.1f%%, %1.1f%%\n',sensitivity(1)*100,specificity(1)*100);
fprintf('AR: %1.1f%%, %1.1f%%\n',sensitivity(2)*100,specificity(2)*100);
fprintf('Paralysis: %1.1f%%, %1.1f%%\n',sensitivity(3)*100,specificity(3)*100);
fprintf('Sensitivity p-value: %1.3f\n',pval.sens);
fprintf('Specificity p-value: %1.3f\n',pval.spec);


%% Compare confidence of reads
fprintf('\n\nAverage percent reads that are confident:\n');
fprintf('Baseline: %1.1f%%\n',sum(method(1).conf(:))/n_reads*100);
fprintf('AR: %1.1f%%\n',sum(method(2).conf(:))/n_reads*100);
fprintf('Paralysis: %1.1f%%\n',sum(method(3).conf(:))/n_reads*100);
fprintf('Confidence p-value: %1.3f\n',pval.conf);

%% Get individual reviewer stats
rev = individual_reviewer_stats(all_ep_or_no,method(4).dc);


%% Breakdown by artifact type
emg = method(4).artifact == 1;
emg_no = method(4).artifact == 0;
fprintf('\nThere was lots of EMG in %d and little EMG for %d\n',...
    sum(emg),sum(emg_no));
n_raters = size(method(1).correct,2);


% Accuracies
fprintf('\n\nAverage percent reads that are correct:\n');
fprintf('Baseline: %1.1f%% for lots emg, %1.1f%% for little emg\n',...
    sum(sum(method(1).correct(emg,:)))/n_raters/sum(emg)*100,...
    sum(sum(method(1).correct(emg_no,:)))/n_raters/sum(emg_no)*100);
fprintf('AR: %1.1f%% for lots emg, %1.1f%% for little emg\n',...
    sum(sum(method(2).correct(emg,:)))/n_raters/sum(emg)*100,...
    sum(sum(method(2).correct(emg_no,:)))/n_raters/sum(emg_no)*100);
fprintf('Paralysis: %1.1f%% for lots emg, %1.1f%% for little emg\n',...
    sum(sum(method(3).correct(emg,:)))/n_raters/sum(emg)*100,...
    sum(sum(method(3).correct(emg_no,:)))/n_raters/sum(emg_no)*100);


% Accuracy increase
fprintf('\n\nIncrease in accuracy from baseline to AR:\n');
fprintf('Lots EMG: %1.1f%%, Little EMG: %1.1f%%\n',...
    (sum(sum(all_acc(:,emg,2),2),1)-sum(sum(all_acc(:,emg,1),2),1))/n_raters/sum(emg)*100,...
    (sum(sum(all_acc(:,emg_no,2),2),1)-sum(sum(all_acc(:,emg_no,1),2),1))/n_raters/sum(emg_no)*100);

% two sample bootstrap
p_val_two = two_sample_bootstrap(all_ep_or_no,method(4).dc,method(4).artifact,all_acc,1e4);
fprintf('P value comparing AR accuracy increase for lots to little EMG: %1.3f\n',p_val_two.acc_ar);
%fprintf('P value comparing AR kappa for twitch to EMG: %1.3f\n',p_val_two.kappa_ar);

%% Make some plots
% Plot % correct for each reviewer by method
%plot_correct(all_acc,results_folder,pval.acc)
%emg_effect_plot(all_acc,emg,emg_no,p_val_two.acc_ar,results_folder)

