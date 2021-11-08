function kappa = fleiss_kappa(reads)

%% Calculate Fleiss Kappa
n_eegs = size(reads,2);
n_reviewers = size(reads,1);
kappa = zeros(3,1);

% Loop over all methods
for m = 1:3

    % Current method
    curr_method = reads(:,:,m); % n_reviewers  x n_eegs
    
    
    % Calculate proportion of all assignments in each category (1 = sz; 0 =  no
    % sz)
    p_sz = sum(sum(curr_method==1))/(n_eegs*n_reviewers);
    p_no = sum(sum(curr_method==0))/(n_eegs*n_reviewers);
    
    % Calculate how much raters agree for each eeg
    p_eegs = zeros(n_eegs,1);
    
    for e = 1:n_eegs % Loop over eegs
        temp_p = 0;
        for k = 0:1 % Loop over no sz or sz
            num_raters_assigning_to_cat = sum(curr_method(:,e) == k);
            
            % increase index of agreement
            temp_p = temp_p + num_raters_assigning_to_cat*(num_raters_assigning_to_cat-1);
        end
        temp_p = 1/(n_reviewers*(n_reviewers-1))*(temp_p);
        p_eegs(e) = temp_p;
    end
    
    % Get mean of the p_eegs
    P_mean = 1/n_eegs*sum(p_eegs);
    
    % Get Pe, chance agreement
    P_e = p_sz^2 + p_no^2;
    
    kappa(m) = (P_mean - P_e)/(1-P_e);

end


end