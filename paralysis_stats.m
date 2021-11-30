clear

%{
I added some papers to paperpile suggesting that a single stage bootstrap
analysis is a reasonable way to handle clustered data. The main paper is 
Field, Christopher A., and Alan H. Welsh. "Bootstrapping clustered data." 
Journal of the Royal Statistical Society: Series B (Statistical Methodology) 
69.3 (2007): 369-390.

This says that the cluster bootstrap method gives  the cluster bootstrap 
gives consistent variance estimates under both the transformation and the 
random-effect model. In the cluster bootstrap method, clusters
are selected by simple random sampling with replacement and there is no 
subsequent permutation.

Another paper, Vanbelle, Sophie, and Adelin Albert. "A bootstrap method 
for comparing correlated kappa coefficients." Journal of Statistical 
Computation and Simulation 78.11 (2008): 1009-1015. This is the reference
justifying the choice to define the one-sided p-value to be the proportion 
of statistic differences less than or equal to zero.


%}

addpath(genpath('.'))

%% Parameters
two_stage = 0; % I think this should be zero
data_folder = '../data/';
file_name = 'Persyst Data Raw 11-30-21 Erin edits (post artifact adjudication)';%'Persyst Data Raw 11-7-21 Erin edits';%';'Persyst Data Raw 02-04-21.xls';
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
mean_stat = [sum(method(1).correct(:))/n_reads*100;...
    sum(method(2).correct(:))/n_reads*100;...
    sum(method(3).correct(:))/n_reads*100];
fprintf('\n\nAverage percent reads that are correct:\n');
write_stats(mean_stat,pval.acc,1)


%% Compare kappa statistics
kappa = fleiss_kappa(all_ep_or_no);
% I verified my algorithm of getting fleiss kappa by building another
% version (alt_fleiss_kappa) that first builds a kappa table in the same
% form as the example table on wikipedia. I tested it by putting in the
% wikipedia table and it gave the same output

main_kappa = kappa;
fprintf('\n\nKappa statistics:\n');
write_stats(kappa,pval.kappa,2)


%% Compare sensitivity and specificity
[sensitivity,specificity] = sens_and_spec(all_ep_or_no,method(4).dc);
fprintf('\n\nSensitivity:\n');
write_stats(sensitivity,pval.sens,2)

fprintf('\n\nSpecificity:\n');
write_stats(specificity,pval.spec,2)

%% Compare confidence of reads
fprintf('\n\nAverage percent reads that are confident:\n');
mean_stat = [sum(method(1).conf(:))/n_reads*100;...
    sum(method(2).conf(:))/n_reads*100;...
    sum(method(3).conf(:))/n_reads*100];
write_stats(mean_stat,pval.conf,1)

%% Get individual reviewer stats
rev = individual_reviewer_stats(all_ep_or_no,method(4).dc);

%% Breakdown by artifact type
emg = method(4).artifact == 1;
emg_no = method(4).artifact == 0;
fprintf('\n\nThere was lots of EMG in %d and little EMG for %d\n',...
    sum(emg),sum(emg_no));
n_raters = size(method(1).correct,2);

%% Post-hoc analysis, is the increase in accuracy for EMG-heavy EEGs significant?
fprintf('\n\nResults for post-hoc analysis looking at EMG-heavy EEGs (%d of %d):\n',...
    sum(emg),length(emg));
% Get the patients with lots of emg
sub_all_ep_or_no = all_ep_or_no(:,emg,:);
sub_dc = method(4).dc(emg);
sub_all_conf = all_conf(:,emg,:);
n_reads_sub = sum(emg)*size(method(1).correct,2);

% Do the bootstrap analysis
emg_heavy_pval = bootstrap_stats(sub_all_ep_or_no,sub_dc,sub_all_conf,1e4,two_stage);

% Compare percentages correct for post-hoc EMG-heavy analysis
mean_stat = [sum(sum(method(1).correct(emg,:)))/n_reads_sub*100;...
    sum(sum(method(2).correct(emg,:)))/n_reads_sub*100;...
    sum(sum(method(3).correct(emg,:)))/n_reads_sub*100];
fprintf('\n\nAverage percent reads that are correct for EMG-heavy EEGs:\n');
write_stats(mean_stat,emg_heavy_pval.acc,1)

% Compare percentages correct for post-hoc EMG-heavy analysis
emg_heavy_kappa = fleiss_kappa(sub_all_ep_or_no);
fprintf('\n\nKappa statistics for EMG-heavy EEGs:\n');
write_stats(emg_heavy_kappa,emg_heavy_pval.kappa,2);

% Compare sensitivity and specificity for post-hoc EMG-heavy analysis
[emg_heavy_sensitivity,emg_heavy_specificity] = sens_and_spec(sub_all_ep_or_no,sub_dc);
fprintf('\n\nSensitivity for EMG-heavy EEGs:\n');
write_stats(emg_heavy_sensitivity,emg_heavy_pval.sens,2)

fprintf('\n\nSpecificity for EMG-heavy:\n');
write_stats(emg_heavy_specificity,emg_heavy_pval.spec,2)

% Compare confidence of reads for post-hoc EMG heavy analysis
fprintf('\n\nAverage percent reads that are confident for EMG-heavy EEGs:\n');
mean_stat = [sum(sum(method(1).conf(emg,:)))/n_reads_sub*100;...
    sum(sum(method(2).conf(emg,:)))/n_reads_sub*100;...
    sum(sum(method(3).conf(emg,:)))/n_reads_sub*100];
write_stats(mean_stat,emg_heavy_pval.conf,1)

%% Post-hoc analysis, is the increase in accuracy for EMG-lite EEGs significant?
fprintf('\n\nResults for post-hoc analysis looking at EMG-lite EEGs (%d of %d):\n',...
    sum(~emg),length(emg));
% Get the patients with little of emg
sub_all_ep_or_no = all_ep_or_no(:,~emg,:);
sub_dc = method(4).dc(~emg);
sub_all_conf = all_conf(:,~emg,:);
n_reads_sub = sum(~emg)*size(method(1).correct,2);

% Do the bootstrap analysis
emg_lite_pval = bootstrap_stats(sub_all_ep_or_no,sub_dc,sub_all_conf,1e4,two_stage);

% Compare percentages correct for post-hoc EMG-lite analysis
mean_stat = [sum(sum(method(1).correct(~emg,:)))/n_reads_sub*100;...
    sum(sum(method(2).correct(~emg,:)))/n_reads_sub*100;...
    sum(sum(method(3).correct(~emg,:)))/n_reads_sub*100];
fprintf('\n\nAverage percent reads that are correct for EMG-lite EEGs:\n');
write_stats(mean_stat,emg_lite_pval.acc,1)

% Compare percentages correct for post-hoc EMG-heavy analysis
emg_lite_kappa = fleiss_kappa(sub_all_ep_or_no);
fprintf('\n\nKappa statistics for EMG-lite EEGs:\n');
write_stats(emg_lite_kappa,emg_lite_pval.kappa,2);

% Compare sensitivity and specificity for post-hoc EMG-heavy analysis
[emg_lite_sensitivity,emg_lite_specificity] = sens_and_spec(sub_all_ep_or_no,sub_dc);
fprintf('\n\nSensitivity for EMG-lite EEGs:\n');
write_stats(emg_lite_sensitivity,emg_lite_pval.sens,2)

fprintf('\n\nSpecificity for EMG-lite:\n');
write_stats(emg_lite_specificity,emg_lite_pval.spec,2)

% Compare confidence of reads for post-hoc EMG heavy analysis
fprintf('\n\nAverage percent reads that are confident for EMG-lite EEGs:\n');
mean_stat = [sum(sum(method(1).conf(~emg,:)))/n_reads_sub*100;...
    sum(sum(method(2).conf(~emg,:)))/n_reads_sub*100;...
    sum(sum(method(3).conf(~emg,:)))/n_reads_sub*100];
write_stats(mean_stat,emg_lite_pval.conf,1)

%% Compare increase in stuff between EMG heavy and lite
% Accuracy increase
fprintf('\n\nDoes AR help more for EMG-heavy EEGs relative to EMG-lite EEGs?:\n');
fprintf('\nIncrease in accuracy from baseline to AR:\n');
fprintf('Lots EMG: %1.1f%%, Little EMG: %1.1f%%\n',...
    (sum(sum(all_acc(:,emg,2),2),1)-sum(sum(all_acc(:,emg,1),2),1))/n_raters/sum(emg)*100,...
    (sum(sum(all_acc(:,emg_no,2),2),1)-sum(sum(all_acc(:,emg_no,1),2),1))/n_raters/sum(emg_no)*100);

fprintf('\nIncrease in kappa from baseline to AR:\n');
fprintf('Lots EMG: %1.3f, Little EMG: %1.3f\n',...
    emg_heavy_kappa(2)-emg_heavy_kappa(1),...
    emg_lite_kappa(2)-emg_lite_kappa(1));

% two sample bootstrap
p_val_two = two_sample_bootstrap(all_ep_or_no,method(4).dc,method(4).artifact,all_acc,1e4);
fprintf('P value comparing AR accuracy increase for lots to little EMG: %1.3f\n',p_val_two.acc_ar);
fprintf('P value comparing AR kappa increase for lots to little EMG: %1.3f\n',p_val_two.kappa_ar);




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

%% All EEGs plot
figure
set(gcf,'position',[222 60 600 737])
tiledlayout(2,1,'tilespacing','tight','padding','compact')

% Accuracy
nexttile
ptitle = 'Average accuracy for each reviewer: all EEGs';
pylabel = 'Percent correct reads';
plot_correct(all_acc,main_pval.acc,ptitle,pylabel)

% Kappa
nexttile
ptitle = 'Interrater reliability: all EEGs';
pylabel = 'Fleiss'' Kappa';
plot_kappas(main_kappa,main_pval.kappa,ptitle,pylabel)
print(gcf,[results_folder,'perc_correct_all'],'-dpng')

%% EMG-heavy vs lite EEGs plot
figure
set(gcf,'position',[222 60 900 737])
tiledlayout(2,2,'tilespacing','compact','padding','compact')

% Accuracy EMG heavy
nexttile
ptitle = 'Reviewer accuracy: EMG-heavy EEGs';
pylabel = 'Percent correct reads';
plot_correct(all_acc(:,emg,:),emg_heavy_pval.acc,ptitle,pylabel)

% Accuracy EMG lite
nexttile
ptitle = 'Reviewer accuracy: EMG-lite EEGs';
pylabel = 'Percent correct reads';
plot_correct(all_acc(:,~emg,:),emg_lite_pval.acc,ptitle,pylabel)

% Kappa EMG heavy
nexttile
ptitle = 'Interrater reliability: EMG-heavy EEGs';
pylabel = 'Fleiss'' Kappa';
plot_kappas(emg_heavy_kappa,emg_heavy_pval.kappa,ptitle,pylabel)

% Kappa EMG lite
nexttile
ptitle = 'Interrater reliability: EMG-lite EEGs';
pylabel = 'Fleiss'' Kappa';
plot_kappas(emg_lite_kappa,emg_lite_pval.kappa,ptitle,pylabel)

print(gcf,[results_folder,'perc_correct_emg'],'-dpng')

%emg_effect_plot(all_acc,emg,emg_no,p_val_two.acc_ar,results_folder)

