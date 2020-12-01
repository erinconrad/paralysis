function paralysis_stats

%% Parameters
data_folder = '../data/';
file_name = 'Persyst Data Raw 11-30-20.xls';
results_folder = '../results/';

%% Get data
method = load_paralysis_data(data_folder,file_name);

%% Convert epileptiform/no-epileptiform to correct/incorrect
for i = 1:3 % Loop over methods
    method(i).correct = is_read_correct(method(i).dc,method(4).dc); % method(4) is report
end

%% Concatenate into one long vector (treating all reads as independent)
for i = 1:3 % Loop over methods
    method(i).correct_vec = method(i).correct(:);
    method(i).conf_vec = method(i).conf(:);
end



%% Convert into one matrix
% Convert into one single 3 dimensional matrix n_reviewers  x n_eegs x
% n_methods
all_reads = nan(size(method(1).correct,2),size(method(1).correct,1),3);
for i = 1:3 % loop over methods
    all_reads(:,:,i) = method(i).correct'; % transpose to convert to n_reviewers x n_eegs
end
n_eegs = size(all_reads,2);
n_reviewers = size(all_reads,1);

%% Basic info

% How many patients are seizing?
fprintf('\n\n%d of %d patients had epileptiform discharges per clinical report.\n',...
    sum(method(4).dc),n_eegs);


%% Compare percentages correct
n_reads = length(method(1).correct(:));
fprintf('\n\nAverage percent reads that are correct:\n');
fprintf('Baseline: %1.1f%%\n',sum(method(1).correct(:))/n_reads*100);
fprintf('AR: %1.1f%%\n',sum(method(2).correct(:))/n_reads*100);
fprintf('Paralysis: %1.1f%%\n',sum(method(3).correct(:))/n_reads*100);

%% Are reads more often correct for each compared to baseline by McNemar test
p12 = mcnemar_test([method(1).correct_vec,method(2).correct_vec]);
p13 = mcnemar_test([method(1).correct_vec,method(3).correct_vec]);
fprintf('\n\nMcNemar test comparing baseline to AR: p = %1.3f\n',p12);
fprintf('McNemar test comparing baseline to paralysis: p = %1.3f\n',p13);

%% Compare kappa statistics
kappa = fleiss_kappa(all_reads);
fprintf('\n\nKappa statistics:\n');
fprintf('Baseline: %1.1f\n',kappa(1));
fprintf('AR: %1.1f\n',kappa(2));
fprintf('Paralysis: %1.1f\n',kappa(3));

%% Compare confidence of reads
fprintf('\n\nAverage percent reads that are confident:\n');
fprintf('Baseline: %1.1f%%\n',sum(method(1).conf(:))/n_reads*100);
fprintf('AR: %1.1f%%\n',sum(method(2).conf(:))/n_reads*100);
fprintf('Paralysis: %1.1f%%\n',sum(method(3).conf(:))/n_reads*100);

%% Are reads more often confident for each compared to baseline by McNemar test
p12 = mcnemar_test([method(1).conf_vec,method(2).conf_vec]);
p13 = mcnemar_test([method(1).conf_vec,method(3).conf_vec]);
fprintf('\n\nMcNemar test comparing baseline to AR: p = %1.3f\n',p12);
fprintf('McNemar test comparing baseline to paralysis: p = %1.3f\n',p13);

%% Make some plots
% Plot % correct for each reviewer by method
plot_correct(all_reads,results_folder)

end