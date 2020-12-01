function get_kappa

%% Use fake data
% Parameters for fake data
fake_perc_correct = [0.6 0.7 0.9];
fake_n_eegs = 10;
fake_n_reviewers = 5;
fake_std_dev = 0.3;

num_correct = generate_fake_data(fake_n_reviewers,fake_n_eegs,fake_perc_correct,fake_std_dev,'kappa');
num_correct = logical(num_correct);

%% Calculate Fleiss Kappa
n_eegs = size(num_correct,2);
n_reviewers = size(num_correct,1);
kappa = zeros(3,1);

% Loop over all methods
for m = 1:3

    % Current method
    curr_method = num_correct(:,:,m); % n_reviewers  x n_eegs
    
    
    % Calculate proportion of all assignments in each category (1 = sz; 0 =  no
    % sz)
    p_sz = sum(sum(curr_method))/(n_eegs*n_reviewers);
    p_no = sum(sum(curr_method==0))/(n_eegs*n_reviewers);
    
    % Calculate how much raters agree for each eeg
    p_eegs = zeros(n_eegs,1);
    for e = 1:n_eegs
        temp_p = 0;
        for k = 0:1
            num_raters_assigning_to_cat = sum(curr_method(:,e) == k);
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