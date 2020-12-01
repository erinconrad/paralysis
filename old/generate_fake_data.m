function num_correct = generate_fake_data(n_reviewers,n_eegs,percs_correct,std_dev,which_method)

if strcmp(which_method,'mcnemar')
    
    % Initialize output
    num_correct = zeros(n_reviewers*n_eegs,3);

    % Loop through method (original, AR, paralysis)
    for j = 1:3

        % Loop through reviewer
        for i = 1:n_reviewers*n_eegs

            % Generate a random percent correct
            perc = percs_correct(j) + randn * std_dev;

            % Set to binary 1 or zero depending on if > 0.5
            if perc > 0.5

                num_correct(i,j) = 1;
                
            end

        end
    end
    
elseif strcmp(which_method,'kappa')
    
    % Initialize output
    num_correct = zeros(n_reviewers,n_eegs,3);
    
    % Loop through method (original, AR, paralysis)
    for j = 1:3

        % Loop through reviewer
        for i = 1:n_reviewers
            
            % Loop through eegs
            for k = 1:n_eegs
                % Generate a random percent correct
                perc = percs_correct(j) + randn * std_dev;

                % Set to binary 1 or zero depending on if > 0.5
                if perc > 0.5
                    
                   num_correct(i,k,j) = 1; 
                end
            end
            
        end
        
    end
    
elseif strcmp(which_method,'kappa2')
    
    % Initialize output
    num_correct = zeros(n_reviewers,n_eegs,3);
    
    for j = 1:3
        for k = 1:n_eegs
            coin = randi(2);
            if coin == 2
                num_correct(:,k,j) = 0;
            else
                num_correct(:,k,j) = 1;
            end
            
            % Loop through reviewers and flip some
            for i = 1:n_reviewers
                coin = randi(10);
                if coin > 9
                    if num_correct(i,k,j) == 1
                        num_correct(i,k,j) = 0;
                    else
                        num_correct(i,k,j) = 1;
                    end
                end
            end
        end
    end
    
else
 
    % Initialize output
    num_correct = zeros(n_reviewers,3);

    % Loop through method (original, AR, paralysis)
    for j = 1:3

        % Loop through reviewer
        for i = 1:n_reviewers

            % Generate a random percent correct
            perc = percs_correct(j) + randn * std_dev;

            % Constrain to be between 0 and 1
            perc = max(0,perc);
            perc = min(1,perc);

            % Get number correct
            num_cor = round(n_eegs * perc);

            num_correct(i,j) = num_cor;

        end
    end
end

end
