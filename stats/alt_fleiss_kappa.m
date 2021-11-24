function kappa = alt_fleiss_kappa(reads)

%% Calculate Fleiss Kappa
kappa = zeros(3,1);

% Loop over all methods
for m = 1:3

    % Current method
    curr_method = reads(:,:,m); % n_reviewers  x n_eegs
    
    % Build kappa table
    kappa_table = build_kappa_table(curr_method);

    % Get kappa from the table
    kappa(m) = kappa_from_table(kappa_table);

end


end