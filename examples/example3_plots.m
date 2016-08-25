% Note: Only works in MATLAB 2014a or later since "histogram" function
% was introduced in that version of MATLAB.
se1=matfile('example3_se1.mat');
se2=matfile('example3_se2.mat');

figure(1);
subplot(2,2,1);
scatter(mean(se1.betasam)', mean(se2.betasam)','b','*');
% Avoid using Statistics Toolbox:
hold on
plot([-1 2],[-1 2],'-','Color','red');
hold off
% ozline = refline([1 0]);
% ozline.Color = 'r';
xlabel('using se_1','Interpreter','none');
ylabel('using se_2','Interpreter','none');
title('Posterior mean of \beta');

subplot(2,2,2);
h1=histogram(se1.hsam,'Normalization','pdf'); 
hold on; 
h2= histogram(se2.hsam,'Normalization','pdf');
title('Posterior density of h');

subplot(2,2,3);
h1=histogram(se1.logpisam,'Normalization','pdf'); 
hold on; 
h2= histogram(se2.logpisam,'Normalization','pdf');
title('Posterior density of log(\pi)');

subplot(2,2,4);
h1=histogram(se1.pvesam,'Normalization','pdf'); 
hold on; 
h2=histogram(se2.pvesam,'Normalization','pdf');
title('Posterior density of PVE');