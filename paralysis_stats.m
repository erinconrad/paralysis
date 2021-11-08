clear

%{
I added some papers to paperpile suggesting that a single stage bootstrap
analysis is a reasonable way to handle clustered data.


%}

%% Parameters
two_stage = 0; % I think this should be zero
data_folder = '../data/';
file_name = 'Persyst Data Raw 11-7-21 Erin edits';%';'Persyst Data Raw 02-04-21.xls';
results_folder = '../results/';

%% seed a rng so same result each time
rng(0)

%% Get data
method = load_paralysis_data(data_folder,file_name);

%% Convert epileptiform/no-epileptiform to correct/incorrect
for i = 1:3 % Loop over methods
    method(i).correct = is_read_correct(method(i).dc,method(4).dc); % method(4) is report
end


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
main_pval = pval;

%% Compare percentages correct
n_reads = length(method(1).correct(:));
fprintf('\n\nAverage percent reads that are correct:\n');
fprintf('Baseline: %1.1f%%\n',sum(method(1).correct(:))/n_reads*100);
fprintf('AR: %1.1f%%\n',sum(method(2).correct(:))/n_reads*100);
fprintf('Paralysis: %1.1f%%\n',sum(method(3).correct(:))/n_reads*100);
fprintf('Accuracy p-value baseline vs AR: %1.3f\n',pval.acc(1));
fprintf('Accuracy p-value baseline vs paralysis: %1.3f\n',pval.acc(2));


%% Compare kappa statistics
kappa = fleiss_kappa(all_ep_or_no);
fprintf('\n\nKappa statistics:\n');
fprintf('Baseline: %1.2f\n',kappa(1));
fprintf('AR: %1.2f\n',kappa(2));
fprintf('Paralysis: %1.2f\n',kappa(3));
fprintf('Kappa p-value baseline vs AR: %1.3f\n',pval.kappa(1));
fprintf('Kappa p-value baseline vs paralysis: %1.3f\n',pval.kappa(2));


%% Compare sensitivity and specificity
[sensitivity,specificity] = sens_and_spec(all_ep_or_no,method(4).dc);
fprintf('\n\nSensitivity and specificity (respectively):\n');
fprintf('Baseline: %1.1f%%, %1.1f%%\n',sensitivity(1)*100,specificity(1)*100);
fprintf('AR: %1.1f%%, %1.1f%%\n',sensitivity(2)*100,specificity(2)*100);
fprintf('Paralysis: %1.1f%%, %1.1f%%\n',sensitivity(3)*100,specificity(3)*100);
fprintf('Sensitivity p-value baseline vs AR: %1.3f\n',pval.sens(1));
fprintf('Sensitivity p-value baseline vs paralysis: %1.3f\n\n',pval.sens(2));
fprintf('Specificity p-value baseline vs AR: %1.3f\n',pval.spec(1));
fprintf('Specificity p-value baseline vs paralysis: %1.3f\n',pval.spec(2));


%% Compare confidence of reads
fprintf('\n\nAverage percent reads that are confident:\n');
fprintf('Baseline: %1.1f%%\n',sum(method(1).conf(:))/n_reads*100);
fprintf('AR: %1.1f%%\n',sum(method(2).conf(:))/n_reads*100);
fprintf('Paralysis: %1.1f%%\n',sum(method(3).conf(:))/n_reads*100);
fprintf('Confidence p-value baseline vs AR: %1.3f\n',pval.conf(1));
fprintf('Confidence p-value baseline vs paralysis: %1.3f\n',pval.conf(2));

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

%% Post-hoc analysis, is the increase in accuracy for EMG-heavy EEGs significant?
fprintf('\nResults for post-hoc analysis looking at EMG-heavy EEGs (%d of %d):\n',...
    sum(emg),length(emg));
% Get the patients with lots of emg
sub_all_ep_or_no = all_ep_or_no(:,emg,:);
sub_dc = method(4).dc(emg);
sub_all_conf = all_conf(:,emg,:);
n_reads_sub = sum(emg)*size(method(1).correct,2);

% Do the bootstrap analysis
pval = bootstrap_stats(sub_all_ep_or_no,sub_dc,sub_all_conf,1e4,two_stage);
emg_heavy_pval = pval;

% Compare percentages correct for post-hoc EMG-heavy analysis
n_reads = length(method(1).correct(:));
fprintf('\n\nAverage percent reads that are correct for EMG-heavy EEGs:\n');
fprintf('Baseline: %1.1f%%\n',sum(sum(method(1).correct(emg,:)))/n_reads_sub*100);
fprintf('AR: %1.1f%%\n',sum(sum(method(2).correct(emg,:)))/n_reads_sub*100);
fprintf('Paralysis: %1.1f%%\n',sum(sum(method(3).correct(emg,:)))/n_reads_sub*100);
fprintf('Accuracy p-value baseline vs AR: %1.3f\n',pval.acc(1));
fprintf('Accuracy p-value baseline vs paralysis: %1.3f\n',pval.acc(2));

% Compare percentages correct for post-hoc EMG-heavy analysis
kappa = fleiss_kappa(sub_all_ep_or_no);
fprintf('\n\nKappa statistics for EMG-heavy EEGs:\n');
fprintf('Baseline: %1.2f\n',kappa(1));
fprintf('AR: %1.2f\n',kappa(2));
fprintf('Paralysis: %1.2f\n',kappa(3));
fprintf('Kappa p-value baseline vs AR: %1.3f\n',pval.kappa(1));
fprintf('Kappa p-value baseline vs paralysis: %1.3f\n',pval.kappa(2));


% Compare sensitivity and specificity for post-hoc EMG-heavy analysis
[sensitivity,specificity] = sens_and_spec(sub_all_ep_or_no,sub_dc);
fprintf('\n\nSensitivity and specificity (respectively) for EMG-heavy EEGs:\n');
fprintf('Baseline: %1.1f%%, %1.1f%%\n',sensitivity(1)*100,specificity(1)*100);
fprintf('AR: %1.1f%%, %1.1f%%\n',sensitivity(2)*100,specificity(2)*100);
fprintf('Paralysis: %1.1f%%, %1.1f%%\n',sensitivity(3)*100,specificity(3)*100);
fprintf('Sensitivity p-value baseline vs AR: %1.3f\n',pval.sens(1));
fprintf('Sensitivity p-value baseline vs paralysis: %1.3f\n\n',pval.sens(2));
fprintf('Specificity p-value baseline vs AR: %1.3f\n',pval.spec(1));
fprintf('Specificity p-value baseline vs paralysis: %1.3f\n',pval.spec(2));


% Compare confidence of reads for post-hoc EMG heavy analysis
fprintf('\n\nAverage percent reads that are confident for EMG-heavy EEGs:\n');
fprintf('Baseline: %1.1f%%\n',sum(sum(method(1).conf(emg,:)))/n_reads_sub*100);
fprintf('AR: %1.1f%%\n',sum(sum(method(2).conf(emg,:)))/n_reads_sub*100);
fprintf('Paralysis: %1.1f%%\n',sum(sum(method(3).conf(emg,:)))/n_reads_sub*100);
fprintf('Confidence p-value baseline vs AR: %1.3f\n',pval.conf(1));
fprintf('Confidence p-value baseline vs paralysis: %1.3f\n',pval.conf(2));

%% Does the presence of discharges predict survival to discharge?
survival = method(4).survival;
dc = method(4).dc;

% get 2x2 table
tbl = crosstab(survival,dc);
[~,p,stats] = fishertest(tbl);
or = stats.OddsRatio;

% descriptive stats
n_surv_with_dc = sum(survival & dc);
perc_surv_with_dc = n_surv_with_dc/sum(dc)*100;

n_surv_without_dc = sum(survival & ~dc);
perc_surv_without_dc = n_surv_without_dc/sum(~dc)*100;

fprintf(['\nDoes the presence of discharges (per report) predict '...
    'survival to hospital discharge?:\n'...
'%d of %d (%1.1f%%) patients with discharges survived to leave the hospital.\n'...
'%d of %d (%1.1f%%) patients without discharges survived to leave the hospital\n'...
'OR %1.1f, p = %1.3f.\n'],...
n_surv_with_dc,sum(dc),perc_surv_with_dc,...
n_surv_without_dc,sum(~dc),perc_surv_without_dc,...
or,p);


%% Make some plots
% Plot % correct for each reviewer by method
plot_correct(all_acc,results_folder,main_pval.acc,'all')
plot_correct(all_acc(:,emg,:),results_folder,emg_heavy_pval.acc,'EMG-heavy')
%emg_effect_plot(all_acc,emg,emg_no,p_val_two.acc_ar,results_folder)

