function [pval,chi2] = mcnemar_test(correct)

a = sum(correct(:,1) & correct(:,2));% orig correct, ar correct
b = sum(correct(:,1) & ~correct(:,2));% orig correct, ar incorrect
c = sum(~correct(:,1) & correct(:,2)); % orig incorrect, ar correct
d = sum(~correct(:,1) & ~correct(:,2)); % orig incorrect, ar incorrect

%% Do exact test
%{
I should probably always use this
%}
if c > b
    new_b = c;
    new_c = b;
else
    new_b = b;
    new_c = c;
end

pval = 0;
n = new_b + new_c;

for i = new_b:n
    pval = pval + nchoosek(n,i);
end
pval = 0.5^n*pval*2;

%% do chi2
%{
This is thought to be appropriate only if b + c > 25
%}
chi2 = (b-c)^2/(b+c);
pval2 = 1- chi2cdf(chi2,1);
end