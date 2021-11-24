function plot_kappas(kappas,pval,ptitle,pylabel)

bar(1:3,kappas)
hold on
psig = p_for_sig_bar(pval);
ylim([0 1])
sigstar({[1 2],[1 3]},psig)
yl = ylim;
%ylim([0 yl(2)])
ylim([0 1])
xticks([1:3])
xlim([0 4])
xticklabels({'Baseline','Artifact reduction','Paralyzed'});
ylabel(pylabel);
title(ptitle);
set(gca,'fontsize',20)

end