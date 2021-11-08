function psig = p_for_sig_bar(pval)

psig = nan(length(pval),1);
for i = 1:length(pval)
    p = pval(i);
    if p < 0.001
        psig(i) = 0.001;
    elseif p < 0.01
        psig(i) = 0.01;
    elseif p < 0.05
        psig(i) = 0.05;
    else
        psig(i) = nan;
    end
end

end