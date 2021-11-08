function ptext = pretty_p_text(pval)

ptext = cell(length(pval),1);
for i = 1:length(pval)
    p = pval(i);
    if p < 0.001
        ptext = '***';
    elseif p < 0.01
        ptext = '**';
    elseif p < 0.05
        ptext = '*';
    else
        ptext = 'na';
    end
end

end