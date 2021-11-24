function kappa_table = build_kappa_table(curr_method)

cats = unique(curr_method);
ncats = length(cats);
nreviewers = size(curr_method,1);
neegs = size(curr_method,2);
kappa_table = nan(neegs,ncats);

for ie = 1:neegs
    for ic = 1:ncats
        c = cats(ic);
        kappa_table(ie,ic) = sum(curr_method(:,ie)==c);
    end
end



end