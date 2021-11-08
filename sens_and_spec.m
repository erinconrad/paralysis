function [sensitivity,specificity] = sens_and_spec(all,true)

n_eegs = size(all,2);
n_reviewers = size(all,1);

sensitivity = zeros(3,1);
specificity = zeros(3,1);

%% Make classification tables
for m = 1:3 % loop over methods
    a = 0; % true pos
    b = 0; % false pos
    c = 0; % false neg
    d = 0; % true neg
    for r = 1:n_reviewers % loop over reviewers
        for e = 1:n_eegs % loop over eegs

            % Get element
            element = all(r,e,m);

            if element == 1 && true(e) == 1 % true positive
                a = a + 1;
            elseif element == 1 && true(e) == 0 % false positive
                b = b + 1;
            elseif element == 0 && true(e) == 1 % false negative
                c = c + 1;
            elseif element == 0 && true(e) == 0 % true negative
                d = d + 1;
            end

        end
    end
    
    sensitivity(m) = a/(a+c); %tp/(tp + fn)
    specificity(m) = d/(b+d); %tn/(tn+fp)
    
end
    
end