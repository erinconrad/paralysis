function cmh(all_reads)

% n_reviewers  x n_eegs x n_methods
%{
4 strata - 4 reviewers

%}


num = 0;
denom = 0;

or_num = 0;
or_denom = 0;

n_strata = size(all_reads,1);

% Loop over strata
for i = 1:n_strata
    
    % Get the matrix
    mat = squeeze(all_reads(i,:,1:2)); % 1:2 because just baseline and AR
    
    % Construct my contingency table
    %{
    a = sum(mat(:,1) == 1 & mat(:,2) == 1);% correct for both
    b = sum(mat(:,1) == 1 & mat(:,2) == 0);% correct baseline, incorrect ar
    c = sum(mat(:,1) == 0 & mat(:,2) == 1);% incorrect baseline, correct ar
    d = sum(mat(:,1) == 0 & mat(:,2) == 0);% incorrect both
    %}
    
    % This seems odd because, again, baseline and AR reads not independent,
    % on same EEGs
    
    a = sum(mat(:,1) == 1); % correct for baseline
    b = sum(mat(:,2) == 1); % correct for ar
    c = sum(mat(:,1) == 0); % incorrect for baseline
    d = sum(mat(:,2) == 0); % incorrect for ar
    n1 = a+b; % total number correct
    n2 = c+d; % total number incorrect
    m1 = a+c; % total number baseline
    m2 = b+d; % total number AR
    t = n1+n2;
    if t ~= m1+m2, error('what'); end
    
    tbl = array2table([a,b;c,d],'VariableNames',{'Baseline','AR'},...
        'RowNames',{'Correct','Incorrect'})
    
    num = num + a - n1*m1/t;
    denom = denom + n1*n2*m1*m2/(t^2*(t-1));
    
    or_num = or_num + a*d/t;
    or_denom = or_denom + b*c/t;
    
end

or = or_num/or_denom

num = num^2;
eta = num/denom;
p = 1 - chi2cdf(eta,1)


end