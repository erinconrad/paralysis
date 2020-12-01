function compare_two_tests

%{
This function uses McNemar's test 
https://en.wikipedia.org/wiki/McNemar%27s_test

It treats each EEG read (e.g., EEG from patient #1 read by Erin) as a
separate subject, and then compares tests (original read vs artifact
removal read vs post-paralysis read) to see how each test performs against
the gold standard clinical read
%}

%% Use fake data
% Parameters for fake data
fake_perc_correct = [0.55 0.65 0.9];
fake_n_eegs = 10;
fake_n_reviewers = 5;
fake_std_dev = 0.2;

% Generate random fake data
%{
This is an N X 3 array where N is the number of eegs x the number of
reviewers, and there are 3 columns (one for each test), and each element is
either a 1 (if the EEG read using that test agrees with the gold standard)
or a 0 (if it disagrees).
%}
num_correct = generate_fake_data(fake_n_reviewers,fake_n_eegs,fake_perc_correct,fake_std_dev,'mcnemar');
num_correct = logical(num_correct);

%% Convert to contingency table
%{
First I restrict to just the first 2 tests (ignore post-paralysis)

This is a contingency table of the form [a,b;c,d] in which the a and d
values are the number of EEG reads in which both tests are correct or both
tests are incorrect. b and c are the number of reads for which one test is
correct and the other is incorrect. If b and c are drastically different,
this suggests that one test is better than the other.
%}
a = sum(num_correct(:,1) & num_correct(:,2));% orig correct, ar correct
b = sum(num_correct(:,1) & ~num_correct(:,2));% orig correct, ar incorrect
c = sum(~num_correct(:,1) & num_correct(:,2)); % orig incorrect, ar correct
d = sum(~num_correct(:,1) & ~num_correct(:,2)); % orig incorrect, ar incorrect

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

fprintf('\nChi2 p-value: %1.3f\nExact p-value: %1.3f\n',pval2,pval);

tbl = [a,b;c,d]


end