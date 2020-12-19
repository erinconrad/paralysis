function [p_val,diff_acc,diff_acc_boot] = perm_test(all_reads,do_plot,perm_vecs)
%{
The idea is that I will randomly swap the reads for baseline vs AR, keeping
reviewer constant. I will measure the difference in AR versus baseline
accuracy and compare this to the distribution of differences I get when
swapping baseline and AR

n_reviewers  x n_eegs x n_methods
%}

baseline = all_reads(:,:,1);
ar = all_reads(:,:,2);
n_reads = size(all_reads,1)*size(all_reads,2);

% Get true accuracies
acc_baseline = sum(sum(baseline))/n_reads;
acc_ar = sum(sum(ar))/n_reads;
diff_acc = acc_ar-acc_baseline;

% Get true kappas
kappa_true = fleiss_kappa(all_reads);
kappa_true_diff = kappa_true(2) - kappa_true(1);

nboot = 1e4;
diff_acc_boot = nan(nboot,1);
diff_kappa_boot = nan(nboot,1);

for ib = 1:nboot
    
    fake_baseline = baseline;
    fake_ar = ar;
    
    if perm_vecs == 1
        for j = 1:size(ar,2) % Loop over eegs
            r = randi([0 1]);
            if r == 1 % swap the whole column representing all reviewers for that read
                fake_baseline(:,j) = ar(:,j);
                fake_ar(:,j) = baseline(:,j);
            end
        end
    else
        % Permute entries
        for i = 1:size(ar,1)
            for j = 1:size(ar,2)
                r = randi([0 1]);
                if r == 1 % swap whether they're looking at baseline or AR
                    fake_baseline(i,j) = ar(i,j);
                    fake_ar(i,j) = baseline(i,j); 
                end
            end
        end
    
    end
    
    % Get new accuracies
    acc_fake_baseline = sum(sum(fake_baseline))/n_reads;
    acc_fake_ar = sum(sum(fake_ar))/n_reads;
    
    diff_acc_boot(ib) = acc_fake_ar-acc_fake_baseline;
    
    % Get new kappa
    fake_reads(:,:,1) = fake_baseline;
    fake_reads(:,:,2) = fake_ar;
    fake_reads(:,:,3) = nan;
    kappa = fleiss_kappa(fake_reads);
    kappa_diff = kappa(2) - kappa(1);
    diff_kappa_boot(ib) = kappa_diff;
    
end

% Sort the bootstrap accuracy differences
diff_acc_boot = sort(diff_acc_boot);
n_more_accurate = sum(diff_acc_boot >= diff_acc);
p_val = n_more_accurate/nboot;
if p_val == 0 
    p_val = 1/(nboot+1);
end

if do_plot
    figure
    plot(diff_acc_boot,'o')
    hold on
    plot(get(gca,'xlim'),[diff_acc diff_acc])
end

% Sort the bootstrap kappa differences
diff_kappa_boot = sort(diff_kappa_boot);
n_more_accurate_k = sum(diff_kappa_boot >= kappa_true_diff);
p_val_k = n_more_accurate_k/nboot;
if p_val_k == 0 
    p_val_k = 1/(nboot+1);
end

end