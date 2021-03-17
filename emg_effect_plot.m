function emg_effect_plot(all_acc,emg,emg_no,pval,results_folder)

colors = [...
    0 0.4470 0.7410;...
    0.8500 0.3250 0.0980];

figure
set(gcf,'position',[440 484 751 314])


n_reviewers = size(all_acc,1);
emg_reads = all_acc(:,emg,1:2);
no_emg_reads = all_acc(:,emg_no,1:2);

all_emg = squeeze(mean(emg_reads,[1]))*100;
all_no = squeeze(mean(no_emg_reads,[1]))*100;

emg_mean = mean(all_emg,1);
no_mean = mean(all_no,1);

plot(emg_mean,'color',colors(1,:),'linewidth',2)
hold on
plot(no_mean,'color',colors(2,:),'linestyle','--','linewidth',2);


ea = plot(1+0.1*rand(size(all_emg,1),1),all_emg(:,1),'o','color',colors(1,:),'markersize',10,'linewidth',3);
em = plot([1-0.05 1+0.05],[emg_mean(1) emg_mean(1)],'linewidth',3,'color',colors(1,:));
    
na = plot(1+0.1*rand(size(all_no,1),1),all_no(:,1),'x','color',colors(2,:),'markersize',10,'linewidth',3);
nm = plot([1-0.05 1+0.05],[no_mean(1) no_mean(1)],'linewidth',3,'color',colors(2,:),'linestyle','--');

plot(2+0.1*rand(size(all_emg,1),1),all_emg(:,2),'o','color',colors(1,:),'markersize',10,'linewidth',3)
plot([2-0.05 2+0.05],[emg_mean(2) emg_mean(2)],'linewidth',3,'color',colors(1,:))
    
plot(2+0.1*rand(size(all_no,1),1),all_no(:,2),'x','color',colors(2,:),'markersize',10,'linewidth',3)
plot([2-0.05 2+0.05],[no_mean(2) no_mean(2)],'linewidth',3,'color',colors(2,:),'linestyle','--')

ylim([-15 115])
legend([ea,em,na,nm],{'Substantial EMG EEGs','Substantial EMG mean','Minimal EMG EEGs','Minimal EMG mean'},...
    'location','south','fontsize',20)
ylabel('Percent correct reads');
xticks([1 2])
xticklabels({'Baseline','Artifact reduction'});
set(gca,'fontsize',20)
print(gcf,[results_folder,'emg'],'-depsc')

if 0
mean_emg = squeeze(mean(emg_reads,[1 2]))*100;
mean_no = squeeze(mean(no_emg_reads,[1 2]))*100;
std_emg = squeeze(std(emg_reads,0,[1 2]))*100;
std_no = squeeze(std(no_emg_reads,0,[1 2]))*100;

errorbar(mean_emg,std_emg,'color',colors(1,:),'linewidth',2);
hold on
errorbar(mean_no,std_no,'color',colors(2,:),'linestyle','--','linewidth',2);
xlim([0.5 2.5])
legend('Substantial EMG EEGs','Minimal EMG','fontsize',20,'location','south')
ylabel('Percentage of reads correct');
xticks([1 2])
xticklabels({'Baseline','Artifact reduction'});
set(gca,'fontsize',20)
end

%{
for i = 1:2 % loop over methods
    
    emg_reads = all_acc(:,emg,i);
    no_emg_reads = all_acc(:,emg_no,i);
    
    mean_emg = mean(emg_reads(:));
    mean_no = mean(no_emg_reads(:));
    std_emg = std(emg_reads(:));
    std_no = std(no_emg_reads(:));
    
    errorbar(i-0.1,mean_emg,std_emg,'color',colors(1,:));
    hold on
    
    errorbar(i+0.1,mean_no,std_no,'color',colors(2,:));
end
%}


end