function plot_correct(all_reads,results_folder,pval)


colors = [...
    0 0.4470 0.7410;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880;...
    0.3010 0.7450 0.9330];

marker_types = {'o';'s';'d';'x'};
marker_sizes = [15;15;15;18];

n_eegs = size(all_reads,2);
n_reviewers = size(all_reads,1);

% Make a vector of for the x-position
x_vec = 0:0.1:0.1*(n_reviewers-1);

rev_id = zeros(n_reviewers,1);
rev_nums = cell(n_reviewers,1);
for i = 1:n_reviewers
    rev_nums{i} = sprintf('Reviewer %d',i);
end

figure
set(gcf,'position',[440 484 751 314])
for i = 1:3 % loop over methods
    
    curr_method = all_reads(:,:,i);
    
    % percent correct for each reviewer
    perc_correct = sum(curr_method,2)/n_eegs*100;
    
    
    
    % Plot
    for j = 1:n_reviewers
      rev_id(j) = plot(i+x_vec(j),perc_correct(j),...
            marker_types{j},...
            'MarkerEdgeColor',colors(j,:),...
            'MarkerFaceColor',colors(j,:),...
            'markersize',marker_sizes(j),...
            'linewidth',3);
        hold on
        
        
    end
    
    % Plot the average
    curr_mean = mean(perc_correct);
    pm = plot([i+x_vec(1) i+x_vec(4)],[curr_mean curr_mean],'k','linewidth',3);
   
    
end

ylim([0 100])
xlim([1+x_vec(1)-0.1,3+x_vec(end)+0.1])
xticks([1:3]+repmat(median(x_vec),1,3))
xticklabels({'Baseline','Artifact reduction','Paralyzed'});
ylabel('Percent correct reads');

if 1
if pval < 0.01
    sigstar({[1,2]},[0.01])
elseif pval < 0.05
    sigstar({[1,2]},[0.05])
else
    sigstar({[1,2]},[nan])
end
end
legend([rev_id;pm],[rev_nums;'Mean'],'location','southeast','fontsize',20)

set(gca,'fontsize',20)

print(gcf,[results_folder,'perc_correct'],'-depsc')
end