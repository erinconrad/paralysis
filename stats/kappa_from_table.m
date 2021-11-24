function kappa = kappa_from_table(kappa_table)

% number of eegs
N = size(kappa_table,1);

% number of raters per eeg
n = sum(kappa_table,2);

% these should all be the same (and should be the number of raters)
if any(diff(n)) ~= 0, error('what'); end
n = n(1);

Pi = 1/(n*(n-1))*(sum(kappa_table.*(kappa_table-1),2));
Pbar = 1/N*sum(Pi);

pj = sum(kappa_table,1)/(N*n);
Pe = sum(pj.^2);

kappa = (Pbar-Pe)/(1-Pe);

end