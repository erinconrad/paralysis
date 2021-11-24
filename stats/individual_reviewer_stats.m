function rev = individual_reviewer_stats(all_reads,truth)
n_eegs = size(all_reads,2);
n_reviewers = size(all_reads,1);
n_methods = size(all_reads,3);

% Convert all reads to accuracy


for ir = 1:n_reviewers
    curr_rev = squeeze(all_reads(ir,:,:));
    
    baseline = curr_rev(:,1);
    ar = curr_rev(:,2);
    
    % Get accuracies
    baseline = is_read_correct(baseline,truth);
    ar = is_read_correct(ar,truth);
    
    % Compare baseline to AR with mcnemar test
    [pval,chi2] = mcnemar_test([baseline,ar]);

    % Number of reads switched to correct vs incorrect
    n_switched_correct = sum(baseline == 0 & ar == 1);
    n_switched_incorrect = sum(baseline == 1 & ar == 0);
    
    % Add to struct
    rev(ir).pval = pval;
    rev(ir).switch_corr = n_switched_correct;
    rev(ir).switch_incorr = n_switched_incorrect;
    
end
    

end