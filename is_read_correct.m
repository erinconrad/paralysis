function correct = is_read_correct(reads,true)

%{
Takes individual reviewer reads (epileptiform or no) and turns them into
correct or incorrect based on whether they agree with ground truth.
%}

correct = nan(size(reads));

for i = 1:size(reads,2)
    correct(:,i) = reads(:,i) == true;
end

end