function write_stats(mean_stat,pval,perc_or_no)

if perc_or_no == 1
    fprintf('Baseline: %1.1f%%\n',mean_stat(1));
    fprintf('AR: %1.1f%%\n',mean_stat(2));
    fprintf('Paralysis: %1.1f%%\n',mean_stat(3));
    fprintf('Accuracy p-value baseline vs AR: %1.3f\n',pval(1));
    fprintf('Accuracy p-value baseline vs paralysis: %1.3f\n',pval(2));
elseif perc_or_no == 2
    fprintf('Baseline: %1.3f\n',mean_stat(1));
    fprintf('AR: %1.3f\n',mean_stat(2));
    fprintf('Paralysis: %1.3f\n',mean_stat(3));
    fprintf('Accuracy p-value baseline vs AR: %1.3f\n',pval(1));
    fprintf('Accuracy p-value baseline vs paralysis: %1.3f\n',pval(2));
end

end