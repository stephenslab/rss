% Run this script after running example4.m.
load('example4-theta04-theta2-100trials.mat');

figure(1);
clf;
subplot(2,2,1);
bplot([alpha_n_diff mu_n_diff s_n_diff]);
set(gca,'YLim',[0 0.002]);
set(gca,'XTick',1:3,'XTickLabel',{'\alpha','\mu','s'},'FontSize',14);
ylabel('maximum absolute difference');
title('Null models');

subplot(2,2,2);
bplot([alpha_e_diff mu_e_diff s_e_diff])
set(gca,'YLim',[0 0.002]);
set(gca,'XTick',1:3,'XTickLabel',{'\alpha','\mu','s'},'FontSize',14);
ylabel('maximum absolute difference');
title('Enrichment models');

subplot(2,2,3);
scatter(log10(bf(:, 1)), log10(bf(:, 2)),'b','o');
hold on
plot([-5 30],[-5 30],':','Color','red');
hold off
% ozline = refline([1 0]);
% ozline.Color = 'r';
xlabel('Log10 Bayes factors from VARBVS','Interpreter','none');
ylabel('Log10 Bayes factors from RSS-VARBVSR','Interpreter','none');
title('');
