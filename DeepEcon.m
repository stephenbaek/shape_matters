%% initialize and prepare data
clear all
close all
clc


getd = @(p)path(p, path);
getd('src');
clear getd

prep_data();

caesar.BirthState = reordercats(categorical(caesar.BirthState),{'Midwest', 'Foreign',  'Northeast', 'South', 'West'}); % changing reference group for BirthState to Midwest
caesar.CarMake = reordercats(categorical(caesar.CarMake),{'Economy', 'Luxury'}); % changing reference group for car make to Economy
caesar.CarModel = reordercats(categorical(caesar.CarModel),{'Non-sedan', 'Sedan',});

summary(caesar)


%% Summary statistics (male only) - Table 1
% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
ReportedBMI = caesar.ReportedWeight./(caesar.ReportedHeight*0.001).^2;
AgeSquared = caesar.Age.^2;
Exp = max(caesar.Age - caesar.Education - 6,0);
ExpSquared = Exp.^2;
CarAge = 2001-caesar.CarYear;
WeightError = caesar.ReportedWeight - caesar.Weight;
BMIError = ReportedBMI - BMI;
HeightError = caesar.ReportedHeight - caesar.Stature;


tbl = [caesar, array2table([BMI, ReportedBMI, WeightError, BMIError, HeightError, AgeSquared,  Exp, ExpSquared, CarAge], 'VariableNames', {'BMI', 'ReportedBMI', 'WeightError', 'BMIError', 'HeightError', 'AgeSquared', 'Experience', 'ExperienceSquared', 'CarAge'})];

% creating a subset that satisfies criteria (ex: selecting males, white collar)
 rows = strcmp(tbl.Gender, 'Male');
 tbl = tbl(rows,:);
 
%selecting variables for analysis
   vars = {'FamilyIncome', 'ReportedHeight', 'ReportedWeight',  'Stature', 'Weight', 'BMI', 'ReportedBMI', 'WeightError', 'BMIError', 'HeightError' ...
           'Experience', 'ExperienceSquared', 'Occupation', 'Education', 'NumberOfChildren',...
           'Fitness', 'MaritalStatus', 'Race', 'BirthState',  'Var1', 'Var2' };
 
 tbl = tbl(:, vars);
 
% find rows with missing data
 TF = ismissing(tbl);
 tbl(any(TF,2),:);
 
% remove rows with missing data
 tbl = rmmissing(tbl);

summary(tbl)

close all
figure('pos',[10 10 900 500]);
marker_size = 5;
marker_color = 'k+';
poly_degree = 1;

x = tbl.Var1; y = tbl.Stature; p = polyfit(x, y, poly_degree); f = polyval(p, x); rsq = 1 - sum((y - f).^2)/sum((y - mean(y)).^2); pf = polyval(p, -3:0.1:3);
subplot(2,3,2); scatter(x, y, marker_size, marker_color); set(gca,'TickLabelInterpreter','latex'); title(sprintf('$P_1$ vs. Height ($R^2=%f$)',rsq),'Interpreter','latex'); xlabel('$P_1$','Interpreter','latex'); ylabel('Height (mm)','Interpreter','latex'); axis([-3,3,-Inf,Inf]); hold on; plot(-3:0.1:3, pf, 'linewidth', 2, 'color', [1,0.5,0]);
x = tbl.Var1; y = tbl.Weight; p = polyfit(x, y, poly_degree); f = polyval(p, x); rsq = 1 - sum((y - f).^2)/sum((y - mean(y)).^2); pf = polyval(p, -3:0.1:3);
subplot(2,3,3); scatter(x, y, marker_size, marker_color); set(gca,'TickLabelInterpreter','latex'); title(sprintf('$P_1$ vs. Weight ($R^2=%f$)',rsq),'Interpreter','latex'); xlabel('$P_1$','Interpreter','latex'); ylabel('Weight (kg)','Interpreter','latex'); axis([-3,3,-Inf,Inf]); hold on; plot(-3:0.1:3, pf, 'linewidth', 2, 'color', [1,0.5,0]);
x = tbl.Var1; y = tbl.BMI; p = polyfit(x, y, poly_degree); f = polyval(p, x); rsq = 1 - sum((y - f).^2)/sum((y - mean(y)).^2); pf = polyval(p, -3:0.1:3);
subplot(2,3,1); scatter(x, y, marker_size, marker_color); set(gca,'TickLabelInterpreter','latex'); title(sprintf('$P_1$ vs. BMI ($R^2=%f$)',rsq),'Interpreter','latex'); xlabel('$P_1$','Interpreter','latex'); ylabel('BMI','Interpreter','latex'); axis([-3,3,-Inf,Inf]); hold on; plot(-3:0.1:3, pf, 'linewidth', 2, 'color', [1,0.5,0]);

x = tbl.Var2; y = tbl.Stature; p = polyfit(x, y, poly_degree); f = polyval(p, x); rsq = 1 - sum((y - f).^2)/sum((y - mean(y)).^2); pf = polyval(p, -3:0.1:3);
subplot(2,3,5); scatter(x, y, marker_size, marker_color); set(gca,'TickLabelInterpreter','latex'); title(sprintf('$P_2$ vs. Height ($R^2=%f$)',rsq),'Interpreter','latex'); xlabel('$P_2$','Interpreter','latex'); ylabel('Height (mm)','Interpreter','latex'); axis([-3,3,-Inf,Inf]); hold on; plot(-3:0.1:3, pf, 'linewidth', 2, 'color', [1,0.5,0]);
x = tbl.Var2; y = tbl.Weight; p = polyfit(x, y, poly_degree); f = polyval(p, x); rsq = 1 - sum((y - f).^2)/sum((y - mean(y)).^2); pf = polyval(p, -3:0.1:3);
subplot(2,3,6); scatter(x, y, marker_size, marker_color); set(gca,'TickLabelInterpreter','latex'); title(sprintf('$P_2$ vs. Weight ($R^2=%f$)',rsq),'Interpreter','latex'); xlabel('$P_2$','Interpreter','latex'); ylabel('Weight (kg)','Interpreter','latex'); axis([-3,3,-Inf,Inf]); hold on; plot(-3:0.1:3, pf, 'linewidth', 2, 'color', [1,0.5,0]);
x = tbl.Var2; y = tbl.BMI; p = polyfit(x, y, poly_degree); f = polyval(p, x); rsq = 1 - sum((y - f).^2)/sum((y - mean(y)).^2); pf = polyval(p, -3:0.1:3);
subplot(2,3,4); scatter(x, y, marker_size, marker_color); set(gca,'TickLabelInterpreter','latex'); title(sprintf('$P_2$ vs. BMI ($R^2=%f$)',rsq),'Interpreter','latex'); xlabel('$P_2$','Interpreter','latex'); ylabel('BMI','Interpreter','latex'); axis([-3,3,-Inf,Inf]); hold on; plot(-3:0.1:3, pf, 'linewidth', 2, 'color', [1,0.5,0]);

%% Summary statistics (Female only) - Table 2
% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
ReportedBMI = caesar.ReportedWeight./(caesar.ReportedHeight*0.001).^2;
AgeSquared = caesar.Age.^2;
Exp = max(caesar.Age - caesar.Education - 6,0);
ExpSquared = Exp.^2;
CarAge = 2001-caesar.CarYear;
WeightError = caesar.ReportedWeight - caesar.Weight;
BMIError = ReportedBMI - BMI;
HeightError = caesar.ReportedHeight - caesar.Stature;

tbl = [caesar, array2table([BMI, ReportedBMI, WeightError, BMIError, HeightError, AgeSquared,  Exp, ExpSquared, CarAge], 'VariableNames', {'BMI', 'ReportedBMI', 'WeightError', 'BMIError', 'HeightError', 'AgeSquared', 'Experience', 'ExperienceSquared', 'CarAge'})];

% creating a subset that satisfies criteria (ex: selecting males, white collar)
 rows = strcmp(tbl.Gender, 'Female');
 tbl = tbl(rows,:);

   vars = {'FamilyIncome', 'ReportedHeight', 'ReportedWeight',  'Stature', 'Weight', 'BMI', 'ReportedBMI', 'WeightError', 'BMIError', 'HeightError' ...
           'Experience', 'ExperienceSquared', 'Occupation', 'Education', 'NumberOfChildren',...
           'Fitness', 'MaritalStatus', 'Race', 'BirthState', 'Var1', 'Var2', 'Var3', 'HipCircumference_Maximum', 'WaistCircumference_Pref', 'CupSize', 'BraSize' };
 
 tbl = tbl(:, vars);
 
% find rows with missing data
 TF = ismissing(tbl);
 tbl(any(TF,2),:);
 
% remove rows with missing data
 tbl = rmmissing(tbl);

summary(tbl)

close all
figure('pos',[10 10 900 500]);
marker_size = 5;
marker_color = 'k+';
poly_degree = 1;

x = tbl.Var1; y = tbl.Stature; p = polyfit(x, y, poly_degree); f = polyval(p, x); rsq = 1 - sum((y - f).^2)/sum((y - mean(y)).^2); pf = polyval(p, -3:0.1:3);
subplot(2,3,2); scatter(x, y, marker_size, marker_color); set(gca,'TickLabelInterpreter','latex'); title(sprintf('$P_1$ vs. Height ($R^2=%f$)',rsq),'Interpreter','latex'); xlabel('$P_1$','Interpreter','latex'); ylabel('Height (mm)','Interpreter','latex'); axis([-3,3,-Inf,Inf]); hold on; plot(-3:0.1:3, pf, 'linewidth', 2, 'color', [1,0.5,0]);
x = tbl.Var1; y = tbl.Weight; p = polyfit(x, y, poly_degree); f = polyval(p, x); rsq = 1 - sum((y - f).^2)/sum((y - mean(y)).^2); pf = polyval(p, -3:0.1:3);
subplot(2,3,3); scatter(x, y, marker_size, marker_color); set(gca,'TickLabelInterpreter','latex'); title(sprintf('$P_1$ vs. Weight ($R^2=%f$)',rsq),'Interpreter','latex'); xlabel('$P_1$','Interpreter','latex'); ylabel('Weight (kg)','Interpreter','latex'); axis([-3,3,-Inf,Inf]); hold on; plot(-3:0.1:3, pf, 'linewidth', 2, 'color', [1,0.5,0]);
x = tbl.Var1; y = tbl.BMI; p = polyfit(x, y, poly_degree); f = polyval(p, x); rsq = 1 - sum((y - f).^2)/sum((y - mean(y)).^2); pf = polyval(p, -3:0.1:3);
subplot(2,3,1); scatter(x, y, marker_size, marker_color); set(gca,'TickLabelInterpreter','latex'); title(sprintf('$P_1$ vs. BMI ($R^2=%f$)',rsq),'Interpreter','latex'); xlabel('$P_1$','Interpreter','latex'); ylabel('BMI','Interpreter','latex'); axis([-3,3,-Inf,Inf]); hold on; plot(-3:0.1:3, pf, 'linewidth', 2, 'color', [1,0.5,0]);

x = tbl.Var2; y = tbl.Stature; p = polyfit(x, y, poly_degree); f = polyval(p, x); rsq = 1 - sum((y - f).^2)/sum((y - mean(y)).^2); pf = polyval(p, -3:0.1:3);
subplot(2,3,5); scatter(x, y, marker_size, marker_color); set(gca,'TickLabelInterpreter','latex'); title(sprintf('$P_2$ vs. Height ($R^2=%f$)',rsq),'Interpreter','latex'); xlabel('$P_2$','Interpreter','latex'); ylabel('Height (mm)','Interpreter','latex'); axis([-3,3,-Inf,Inf]); hold on; plot(-3:0.1:3, pf, 'linewidth', 2, 'color', [1,0.5,0]);
x = tbl.Var2; y = tbl.Weight; p = polyfit(x, y, poly_degree); f = polyval(p, x); rsq = 1 - sum((y - f).^2)/sum((y - mean(y)).^2); pf = polyval(p, -3:0.1:3);
subplot(2,3,6); scatter(x, y, marker_size, marker_color); set(gca,'TickLabelInterpreter','latex'); title(sprintf('$P_2$ vs. Weight ($R^2=%f$)',rsq),'Interpreter','latex'); xlabel('$P_2$','Interpreter','latex'); ylabel('Weight (kg)','Interpreter','latex'); axis([-3,3,-Inf,Inf]); hold on; plot(-3:0.1:3, pf, 'linewidth', 2, 'color', [1,0.5,0]);
x = tbl.Var2; y = tbl.BMI; p = polyfit(x, y, poly_degree); f = polyval(p, x); rsq = 1 - sum((y - f).^2)/sum((y - mean(y)).^2); pf = polyval(p, -3:0.1:3);
subplot(2,3,4); scatter(x, y, marker_size, marker_color); set(gca,'TickLabelInterpreter','latex'); title(sprintf('$P_2$ vs. BMI ($R^2=%f$)',rsq),'Interpreter','latex'); xlabel('$P_2$','Interpreter','latex'); ylabel('BMI','Interpreter','latex'); axis([-3,3,-Inf,Inf]); hold on; plot(-3:0.1:3, pf, 'linewidth', 2, 'color', [1,0.5,0]);

%% Nayadara-Watson kernel-regression (height, male only) - Figure 4
tbl = caesar;

% creating a subset that satisfies criteria (ex: selecting males, white collar)
rows = strcmp(tbl.Gender, 'Male');
tbl = tbl(rows,:);
 
vars = {'ReportedHeight', 'Stature'};
tbl = tbl(:, vars);
 
% find rows with missing data
TF = ismissing(tbl);
tbl(any(TF,2),:);
 
% remove rows with missing data
tbl = rmmissing(tbl);

y = tbl.ReportedHeight - tbl.Stature;
x = tbl.Stature;

bandwidth = 1.06*std(x)*(size(tbl,1)^(-0.2));
range = 1550:10:2000;

estimator = [];
for x_chosen = range
    num = sum(y.*EpanechnikovKernel( (x - x_chosen)/bandwidth ));
    den = sum(EpanechnikovKernel( (x - x_chosen)/bandwidth ));
    estimator = [estimator, num/den];
end

% bootstrapping
n = size(tbl, 1);
bb = 500;
idx = randi(n, [bb, n]);

x_temp = x(idx);
y_temp = y(idx);

estimator_boot = [];
for x_chosen = range
    num = sum(y_temp.*EpanechnikovKernel( (x_temp - x_chosen)/bandwidth ), 2);
    den = sum(EpanechnikovKernel( (x_temp - x_chosen)/bandwidth ), 2);
    estimator_boot = [estimator_boot, num./den];
end

estimator_boot(find(any(isnan(estimator_boot),2)), :) = [];
estimator_boot = estimator_boot(randperm(250),:);
conf_upper = estimator + 1.96*std(estimator_boot);
conf_lower = estimator - 1.96*std(estimator_boot);

figure('pos',[10 10 720 270]); subplot(1,2,1); hp = plot(range, estimator, 'linewidth', 2, 'color', [0,0,0]); hold on; hc = plot(range, zeros(size(range)), '-k');
hd = plot(range, conf_upper, '--k');
plot(range, conf_lower, '--k');
title('Reporting Error in Height (Male)','Interpreter','latex');
xlabel('Measured Height (mm)','Interpreter','latex');
ylabel('Mean of Reporting Error (mm)','Interpreter','latex');
lg = legend([hp,hd], 'Reporting Error','95\% Confidence Band', 'Location', 'northeast');
set(lg, 'interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex');
axis([range(1), range(end), -20, 70])


%% Nayadara-Watson kernel-regression (height, female only) - Figure 4
tbl = caesar;

% creating a subset that satisfies criteria (ex: selecting males, white collar)
rows = strcmp(tbl.Gender, 'Female');
tbl = tbl(rows,:);
 
vars = {'ReportedHeight', 'Stature'};
tbl = tbl(:, vars);
 
% find rows with missing data
TF = ismissing(tbl);
tbl(any(TF,2),:);
 
% remove rows with missing data
tbl = rmmissing(tbl);

y = tbl.ReportedHeight - tbl.Stature;
x = tbl.Stature;

bandwidth = 1.06*std(x)*(size(tbl,1)^(-0.2));
range = 1450:10:1850;

estimator = [];
for x_chosen = range
    num = sum(y.*EpanechnikovKernel( (x - x_chosen)/bandwidth ));
    den = sum(EpanechnikovKernel( (x - x_chosen)/bandwidth ));
    estimator = [estimator, num/den];
end

% bootstrapping
n = size(tbl, 1);
bb = 500;
idx = randi(n, [bb, n]);

x_temp = x(idx);
y_temp = y(idx);

estimator_boot = [];
for x_chosen = range
    num = sum(y_temp.*EpanechnikovKernel( (x_temp - x_chosen)/bandwidth ), 2);
    den = sum(EpanechnikovKernel( (x_temp - x_chosen)/bandwidth ), 2);
    estimator_boot = [estimator_boot, num./den];
end

estimator_boot(find(any(isnan(estimator_boot),2)), :) = [];
estimator_boot = estimator_boot(randperm(250),:);
conf_upper = estimator + 1.96*std(estimator_boot);
conf_lower = estimator - 1.96*std(estimator_boot);

subplot(1,2,2); hp = plot(range, estimator, 'linewidth', 2, 'color', [0,0,0]); hold on; hc = plot(range, zeros(size(range)), '-k');
hd = plot(range, conf_upper, '--k');
plot(range, conf_lower, '--k');
title('Reporting Error in Height (Female)','Interpreter','latex');
xlabel('Measured Height (mm)','Interpreter','latex');
ylabel('Mean of Reporting Error (mm)','Interpreter','latex');
lg = legend([hp,hd], 'Reporting Error','95\% Confidence Band', 'Location', 'northeast');
set(lg, 'interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex');
axis([range(1), range(end), -20, 70])


%% Nayadara-Watson kernel-regression (weight, male only) - Figure 5
tbl = caesar;

% creating a subset that satisfies criteria (ex: selecting males, white collar)
rows = strcmp(tbl.Gender, 'Male');
tbl = tbl(rows,:);
 
vars = {'ReportedWeight', 'Weight'};
tbl = tbl(:, vars);
 
% find rows with missing data
TF = ismissing(tbl);
tbl(any(TF,2),:);
 
% remove rows with missing data
tbl = rmmissing(tbl);

y = tbl.ReportedWeight - tbl.Weight;
x = tbl.Weight;

bandwidth = 1.06*std(x)*(size(tbl,1)^(-0.2));
range = 45:1:120;

estimator = [];
for x_chosen = range
    num = sum(y.*EpanechnikovKernel( (x - x_chosen)/bandwidth ));
    den = sum(EpanechnikovKernel( (x - x_chosen)/bandwidth ));
    estimator = [estimator, num/den];
end

% bootstrapping
n = size(tbl, 1);
bb = 500;
idx = randi(n, [bb, n]);

x_temp = x(idx);
y_temp = y(idx);

estimator_boot = [];
for x_chosen = range
    num = sum(y_temp.*EpanechnikovKernel( (x_temp - x_chosen)/bandwidth ), 2);
    den = sum(EpanechnikovKernel( (x_temp - x_chosen)/bandwidth ), 2);
    estimator_boot = [estimator_boot, num./den];
end

estimator_boot(find(any(isnan(estimator_boot),2)), :) = [];
estimator_boot = estimator_boot(randperm(250),:);
conf_upper = estimator + 1.96*std(estimator_boot);
conf_lower = estimator - 1.96*std(estimator_boot);

figure('pos',[10 10 720 270]); subplot(1,2,1); hp = plot(range, estimator, 'linewidth', 2, 'color', [0,0,0]); hold on; hc = plot(range, zeros(size(range)), '-k');
hd = plot(range, conf_upper, '--k');
plot(range, conf_lower, '--k');
title('Reporting Error in Weight (Male)','Interpreter','latex');
xlabel('Measured Weight (kg)','Interpreter','latex');
ylabel('Mean of Reporting Error (kg)','Interpreter','latex');
lg = legend([hp,hd], 'Reporting Error','95\% Confidence Band', 'Location', 'northeast');
set(lg, 'interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex');
axis([range(1), range(end), -5.5, 5])

%% Nayadara-Watson kernel-regression (weight, female only) - Figure 5
tbl = caesar;

% creating a subset that satisfies criteria (ex: selecting males, white collar)
rows = strcmp(tbl.Gender, 'Female');
tbl = tbl(rows,:);
 
vars = {'ReportedWeight', 'Weight'};
tbl = tbl(:, vars);
 
% find rows with missing data
TF = ismissing(tbl);
tbl(any(TF,2),:);
 
% remove rows with missing data
tbl = rmmissing(tbl);

y = tbl.ReportedWeight - tbl.Weight;
x = tbl.Weight;

bandwidth = 1.06*std(x)*(size(tbl,1)^(-0.2));
range = 40:1:100;

estimator = [];
for x_chosen = range
    num = sum(y.*EpanechnikovKernel( (x - x_chosen)/bandwidth ));
    den = sum(EpanechnikovKernel( (x - x_chosen)/bandwidth ));
    estimator = [estimator, num/den];
end

% bootstrapping
n = size(tbl, 1);
bb = 500;
idx = randi(n, [bb, n]);

x_temp = x(idx);
y_temp = y(idx);

estimator_boot = [];
for x_chosen = range
    num = sum(y_temp.*EpanechnikovKernel( (x_temp - x_chosen)/bandwidth ), 2);
    den = sum(EpanechnikovKernel( (x_temp - x_chosen)/bandwidth ), 2);
    estimator_boot = [estimator_boot, num./den];
end

estimator_boot(find(any(isnan(estimator_boot),2)), :) = [];
estimator_boot = estimator_boot(randperm(250),:);
conf_upper = estimator + 1.96*std(estimator_boot);
conf_lower = estimator - 1.96*std(estimator_boot);

subplot(1,2,2); hp = plot(range, estimator, 'linewidth', 2, 'color', [0,0,0]); hold on; hc = plot(range, zeros(size(range)), '-k');
hd = plot(range, conf_upper, '--k');
plot(range, conf_lower, '--k');
title('Reporting Error in Weight (Female)','Interpreter','latex');
xlabel('Measured Weight (kg)','Interpreter','latex');
ylabel('Mean of Reporting Error (kg)','Interpreter','latex');
lg = legend([hp,hd], 'Reporting Error','95\% Confidence Band', 'Location', 'northeast');
set(lg, 'interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex');
axis([range(1), range(end), -5.5, 5])



%% Quantile Regression Height (male only) - Figure 6
tbl = caesar;

% creating a subset that satisfies criteria (ex: selecting males, white collar)
rows = strcmp(tbl.Gender, 'Male');
tbl = tbl(rows,:);
 
vars = {'ReportedHeight', 'Stature'};
tbl = tbl(:, vars);
 
% find rows with missing data
TF = ismissing(tbl);
tbl(any(TF,2),:);
 
% remove rows with missing data
tbl = rmmissing(tbl);
 
x = tbl.Stature;
y = tbl.ReportedHeight-tbl.Stature;
[x, I] = sort(x);
y = y(I);

order = 1; 
[p,stats] = quantreg(x, y, 0.1, order);

lower_decile = polyval(p,x);
lower_decile_stats = stats;
[p,stats] = quantreg(x, y, 0.25, order);
lower_quartile = polyval(p,x);
lower_quartile_stats = stats;
[p,stats] = quantreg(x, y, 0.5, order);
median = polyval(p,x);
median_p = p;
median_stats = stats;
[p,stats] = quantreg(x, y, 0.75, order);
upper_quartile = polyval(p,x);
upper_quartile_stats = stats;
[p,stats] = quantreg(x, y, 0.9, order);
upper_decile = polyval(p,x);
upper_decile_stats = stats;

figure('pos',[10 10 720 270]);

subplot(1,2,1); 
hold on
h = patch([x;x(end:-1:1);x(1)],...
    [median_stats.yfitci(:,1);median_stats.yfitci(end:-1:1,2);median_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

h = patch([x;x(end:-1:1);x(1)],...
    [upper_quartile_stats.yfitci(:,1);upper_quartile_stats.yfitci(end:-1:1,2);upper_quartile_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

h = patch([x;x(end:-1:1);x(1)],...
    [lower_quartile_stats.yfitci(:,1);lower_quartile_stats.yfitci(end:-1:1,2);lower_quartile_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

h = patch([x;x(end:-1:1);x(1)],...
    [upper_decile_stats.yfitci(:,1);upper_decile_stats.yfitci(end:-1:1,2);upper_decile_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

h = patch([x;x(end:-1:1);x(1)],...
    [lower_decile_stats.yfitci(:,1);lower_decile_stats.yfitci(end:-1:1,2);lower_decile_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

linewidth = 2;
hd = plot(x,lower_decile, ':k', 'linewidth', linewidth)
hq = plot(x,lower_quartile, '--k', 'linewidth', linewidth)
hm = plot(x,median, '-k', 'linewidth', linewidth)
plot(x,upper_quartile, '--k', 'linewidth', linewidth)
plot(x,upper_decile, ':k', 'linewidth', linewidth)
plot([min(x), max(x)], [0,0], '-k')
axis([1550, 2000, -20, 70]);
lg = legend([hm, hq, hd, h], 'Median', 'Upper/Lower Quartile', 'Upper/Lower Decile', '95\% Confidence Band', 'Location', 'northeast');
set(lg, 'interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex');
xlabel('Measured Height (mm)', 'interpreter', 'latex')
ylabel('Quantiles of Reporting Error (mm)', 'interpreter', 'latex')
title('Reporting Error in Height (Male)', 'interpreter', 'latex')
box on


%% Quantile Regression Height (female only) - Figure 6
tbl = caesar;

% creating a subset that satisfies criteria (ex: selecting males, white collar)
rows = strcmp(tbl.Gender, 'Female');
tbl = tbl(rows,:);
 
vars = {'ReportedHeight', 'Stature'};
tbl = tbl(:, vars);
 
% find rows with missing data
TF = ismissing(tbl);
tbl(any(TF,2),:);
 
% remove rows with missing data
tbl = rmmissing(tbl);
 
x = tbl.Stature;
y = tbl.ReportedHeight - tbl.Stature;
[x, I] = sort(x);
y = y(I);

order = 1;
[p,stats] = quantreg(x, y, 0.1, order);

lower_decile = polyval(p,x);
lower_decile_stats = stats;
[p,stats] = quantreg(x, y, 0.25, order);
lower_quartile = polyval(p,x);
lower_quartile_stats = stats;
[p,stats] = quantreg(x, y, 0.5, order);
median = polyval(p,x);
median_p = p;
median_stats = stats;
[p,stats] = quantreg(x, y, 0.75, order);
upper_quartile = polyval(p,x);
upper_quartile_stats = stats;
[p,stats] = quantreg(x, y, 0.9, order);
upper_decile = polyval(p,x);
upper_decile_stats = stats;

subplot(1,2,2); 
hold on
h = patch([x;x(end:-1:1);x(1)],...
    [median_stats.yfitci(:,1);median_stats.yfitci(end:-1:1,2);median_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

h = patch([x;x(end:-1:1);x(1)],...
    [upper_quartile_stats.yfitci(:,1);upper_quartile_stats.yfitci(end:-1:1,2);upper_quartile_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

h = patch([x;x(end:-1:1);x(1)],...
    [lower_quartile_stats.yfitci(:,1);lower_quartile_stats.yfitci(end:-1:1,2);lower_quartile_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

h = patch([x;x(end:-1:1);x(1)],...
    [upper_decile_stats.yfitci(:,1);upper_decile_stats.yfitci(end:-1:1,2);upper_decile_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

h = patch([x;x(end:-1:1);x(1)],...
    [lower_decile_stats.yfitci(:,1);lower_decile_stats.yfitci(end:-1:1,2);lower_decile_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

linewidth = 2;
hd = plot(x,lower_decile, ':k', 'linewidth', linewidth)
hq = plot(x,lower_quartile, '--k', 'linewidth', linewidth)
hm = plot(x,median, '-k', 'linewidth', linewidth)
plot(x,upper_quartile, '--k', 'linewidth', linewidth)
plot(x,upper_decile, ':k', 'linewidth', linewidth)
plot([min(x), max(x)], [0,0], '-k')
axis([1450, 1850, -20, 70]);
lg = legend([hm, hq, hd, h], 'Median', 'Upper/Lower Quartile', 'Upper/Lower Decile', '95\% Confidence Band', 'Location', 'northeast');
set(lg, 'interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex');
xlabel('Measured Height (mm)', 'interpreter', 'latex')
ylabel('Quantiles of Reporting Error (mm)', 'interpreter', 'latex')
title('Reporting Error in Height (Female)', 'interpreter', 'latex')
box on


%% Quantile Regression Weight (male only) - Figure 7
tbl = caesar;

% creating a subset that satisfies criteria (ex: selecting males, white collar)
rows = strcmp(tbl.Gender, 'Male');
tbl = tbl(rows,:);
 
vars = {'ReportedWeight', 'Weight'};
tbl = tbl(:, vars);
 
% find rows with missing data
TF = ismissing(tbl);
tbl(any(TF,2),:);
 
% remove rows with missing data
tbl = rmmissing(tbl);
 
x = tbl.Weight;
y = tbl.ReportedWeight - tbl.Weight;
[x, I] = sort(x);
y = y(I);

order = 2; 
[p,stats] = quantreg(x, y, 0.1, order);

lower_decile = polyval(p,x);
lower_decile_stats = stats;
[p,stats] = quantreg(x, y, 0.25, order);
lower_quartile = polyval(p,x);
lower_quartile_stats = stats;
[p,stats] = quantreg(x, y, 0.5, order);
median = polyval(p,x);
median_p = p;
median_stats = stats;
[p,stats] = quantreg(x, y, 0.75, order);
upper_quartile = polyval(p,x);
upper_quartile_stats = stats;
[p,stats] = quantreg(x, y, 0.9, order);
upper_decile = polyval(p,x);
upper_decile_stats = stats;

figure('pos',[10 10 720 270]);

subplot(1,2,1); 
hold on
h = patch([x;x(end:-1:1);x(1)],...
    [median_stats.yfitci(:,1);median_stats.yfitci(end:-1:1,2);median_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

h = patch([x;x(end:-1:1);x(1)],...
    [upper_quartile_stats.yfitci(:,1);upper_quartile_stats.yfitci(end:-1:1,2);upper_quartile_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

h = patch([x;x(end:-1:1);x(1)],...
    [lower_quartile_stats.yfitci(:,1);lower_quartile_stats.yfitci(end:-1:1,2);lower_quartile_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

h = patch([x;x(end:-1:1);x(1)],...
    [upper_decile_stats.yfitci(:,1);upper_decile_stats.yfitci(end:-1:1,2);upper_decile_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

h = patch([x;x(end:-1:1);x(1)],...
    [lower_decile_stats.yfitci(:,1);lower_decile_stats.yfitci(end:-1:1,2);lower_decile_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

linewidth = 2;
hd = plot(x,lower_decile, ':k', 'linewidth', linewidth)
hq = plot(x,lower_quartile, '--k', 'linewidth', linewidth)
hm = plot(x,median, '-k', 'linewidth', linewidth)
plot(x,upper_quartile, '--k', 'linewidth', linewidth)
plot(x,upper_decile, ':k', 'linewidth', linewidth)
plot([min(x), max(x)], [0,0], '-k')
axis([55, 120, -10, 5])
lg = legend([hm, hq, hd, h], 'Median', 'Upper/Lower Quartile', 'Upper/Lower Decile', '95\% Confidence Band', 'Location', 'southwest');
set(lg, 'interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex');
xlabel('Measured Weight (kg)', 'interpreter', 'latex')
ylabel('Quantiles of Reporting Error (kg)', 'interpreter', 'latex')
title('Reporting Error in Weight (Male)', 'interpreter', 'latex')
box on

%% Quantile Regression Weight (female only) - Figure 7
tbl = caesar;

% creating a subset that satisfies criteria (ex: selecting males, white collar)
rows = strcmp(tbl.Gender, 'Female');
tbl = tbl(rows,:);
 
vars = {'ReportedWeight', 'Weight'};
tbl = tbl(:, vars);
 
% find rows with missing data
TF = ismissing(tbl);
tbl(any(TF,2),:);
 
% remove rows with missing data
tbl = rmmissing(tbl);
 
x = tbl.Weight;
y = tbl.ReportedWeight - tbl.Weight;
[x, I] = sort(x);
y = y(I);

order = 2; 
[p,stats] = quantreg(x, y, 0.1, order);

lower_decile = polyval(p,x);
lower_decile_stats = stats;
[p,stats] = quantreg(x, y, 0.25, order);
lower_quartile = polyval(p,x);
lower_quartile_stats = stats;
[p,stats] = quantreg(x, y, 0.5, order);
median = polyval(p,x);
median_p = p;
median_stats = stats;
[p,stats] = quantreg(x, y, 0.75, order);
upper_quartile = polyval(p,x);
upper_quartile_stats = stats;
[p,stats] = quantreg(x, y, 0.9, order);
upper_decile = polyval(p,x);
upper_decile_stats = stats;

subplot(1,2,2); 
hold on
h = patch([x;x(end:-1:1);x(1)],...
    [median_stats.yfitci(:,1);median_stats.yfitci(end:-1:1,2);median_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

h = patch([x;x(end:-1:1);x(1)],...
    [upper_quartile_stats.yfitci(:,1);upper_quartile_stats.yfitci(end:-1:1,2);upper_quartile_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

h = patch([x;x(end:-1:1);x(1)],...
    [lower_quartile_stats.yfitci(:,1);lower_quartile_stats.yfitci(end:-1:1,2);lower_quartile_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

h = patch([x;x(end:-1:1);x(1)],...
    [upper_decile_stats.yfitci(:,1);upper_decile_stats.yfitci(end:-1:1,2);upper_decile_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

h = patch([x;x(end:-1:1);x(1)],...
    [lower_decile_stats.yfitci(:,1);lower_decile_stats.yfitci(end:-1:1,2);lower_decile_stats.yfitci(1,1)],...
    'r');
set(h, 'FaceColor', [0.9,0.9,0.9])
set(h, 'EdgeColor', 'None')

linewidth = 2;
hd = plot(x,lower_decile, ':k', 'linewidth', linewidth)
hq = plot(x,lower_quartile, '--k', 'linewidth', linewidth)
hm = plot(x,median, '-k', 'linewidth', linewidth)
plot(x,upper_quartile, '--k', 'linewidth', linewidth)
plot(x,upper_decile, ':k', 'linewidth', linewidth)
plot([min(x), max(x)], [0,0], '-k')
% xlim([40, 100])
axis([40, 100, -10, 5])
lg = legend([hm, hq, hd, h], 'Median', 'Upper/Lower Quartile', 'Upper/Lower Decile', '95\% Confidence Band', 'Location', 'southwest');
set(lg, 'interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex');
xlabel('Measured Weight (kg)', 'interpreter', 'latex')
ylabel('Quantiles of Reporting Error (kg)', 'interpreter', 'latex')
title('Reporting Error in Weight (Female)', 'interpreter', 'latex')
box on



%% LASSO - Various body measures (Male only) - Figure 10
% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
AgeSquared = caesar.Age.^2;
Exp = max(caesar.Age - caesar.Education - 6,0);
ExpSquared = Exp.^2;
CarAge = 2001-caesar.CarYear;

tbl = [caesar, array2table([BMI, AgeSquared, Exp, ExpSquared, CarAge], 'VariableNames', {'BMI', 'AgeSquared','Experience', 'ExperienceSquared', 'CarAge'})];

% creating a subset that satisfies criteria (ex: selecting males, white collar)
gender = 'Male';
rows = strcmp(tbl.Gender, gender);
tbl = tbl(rows,:);
 
% selecting variables for analysis
vars = {'FamilyIncome',  'Experience', 'ExperienceSquared', 'Occupation', 'Education', 'NumberOfChildren',...
        'MaritalStatus', 'Race'};
    
n_others = length(vars);
if strcmp(gender, 'Male') == 1
    vars = cat(2, vars, caesar.Properties.VariableNames(setdiff(26:70,[36,28,29,61,62,63])));
else
    vars = cat(2, vars, caesar.Properties.VariableNames(setdiff(26:70,[28,29,61,62,63])));
end
tbl = tbl(:, vars);
n_measures = length(vars) - n_others;
 

% find rows with missing data
 TF = ismissing(tbl);
 tbl(any(TF,2),:);
 
% remove rows with missing data
 tbl = rmmissing(tbl);
 
% take log on family income
 tbl.FamilyIncome = log(tbl.FamilyIncome);

 % LASSO
tbl2 = tbl;
tbl2 = removevars(tbl, {'FamilyIncome', 'Occupation', 'MaritalStatus', 'Race'});
tbl2 = [array2table([dummyvar(categorical(tbl.Occupation)),...
                           dummyvar(categorical(tbl.MaritalStatus)),...
                           dummyvar(categorical(tbl.Race))],...
         'VariableNames', {'Occupation_Blue', 'Occupation_Mng', 'Occupation_Service', 'Occupation_White',...
                           'MaritalStatus_Divorced_or_Widowed','MaritalStatus_Married','MaritalStatus_Single',...
                           'Race_Asian','Race_Black','Race_Hispanic','Race_White'}),...
         tbl2];
tbl2 = removevars(tbl2, {'Occupation_White', 'MaritalStatus_Single', 'Race_White'});

n_others = size(tbl2,2) - n_measures;
tbl3 = table2array(tbl2(:, n_others+1:end));
tbl3 = (tbl3 - mean(tbl3))./std(tbl3);
D = x2fx(tbl3, 'Interaction');
D(:,1) = [];
D = [table2array(tbl2(:, 1:n_others)), D];

[B, fitInfo] = lasso(D, tbl.FamilyIncome, 'CV', 10);


figure('pos',[10 10 720 270]);

a = subplot(1,2,1);
lassoPlot(B, fitInfo, 'PlotType', 'CV', 'parent', a)
l = legend('show');
set(l, 'Interpreter', 'latex');
set(get(a,'title'), 'Interpreter', 'latex');
set(get(a,'title'), 'String', 'Cross-Validated MSE of Lasso Fit (Male)');
set(get(a,'xlabel'), 'Interpreter', 'latex');
set(get(a,'ylabel'), 'Interpreter', 'latex');
set(a,'TickLabelInterpreter', 'latex');
 
fprintf('==============================================\n');
fprintf('|  LASSO RESULT (MALE)                       |\n');
fprintf('----------------------------------------------\n');
fprintf('|    Lambda @ 1SE: %f\n', fitInfo.Lambda1SE);
fprintf('|    Lambda @ Min: %f\n', fitInfo.LambdaMinMSE);
fprintf('|    SE @ 1SE: %f\n', fitInfo.SE(fitInfo.Index1SE));
fprintf('|    SE @ Min: %f\n', fitInfo.SE(fitInfo.IndexMinMSE));
fprintf('|    MSE @ 1SE: %f\n', fitInfo.MSE(fitInfo.Index1SE));
fprintf('|    MSE @ Min: %f\n', fitInfo.MSE(fitInfo.IndexMinMSE));
fprintf('===============================================\n');

mdl_lasso = fitlm(D(:, find(B(:,fitInfo.IndexMinMSE))), tbl.FamilyIncome);

%% LASSO - Various body measures (Female only) - Figure 10
% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
AgeSquared = caesar.Age.^2;
Exp = max(caesar.Age - caesar.Education - 6,0);
ExpSquared = Exp.^2;
CarAge = 2001-caesar.CarYear;

tbl = [caesar, array2table([BMI, AgeSquared, Exp, ExpSquared, CarAge], 'VariableNames', {'BMI', 'AgeSquared','Experience', 'ExperienceSquared', 'CarAge'})];

% creating a subset that satisfies criteria (ex: selecting males, white collar)
gender = 'Female';
rows = strcmp(tbl.Gender, gender);
tbl = tbl(rows,:);
 
% selecting variables for analysis
vars = {'FamilyIncome',  'Experience', 'ExperienceSquared', 'Occupation', 'Education', 'NumberOfChildren',...
       'Fitness', 'MaritalStatus', 'Race', 'BirthState'};
    
n_others = length(vars);
if strcmp(gender, 'Male') == 1
    vars = cat(2, vars, caesar.Properties.VariableNames(setdiff(26:70,[36,28,29,61,62,63])));
else
    vars = cat(2, vars, caesar.Properties.VariableNames(setdiff(26:70,[28,29,61,62,63])));
end
tbl = tbl(:, vars);
n_measures = length(vars) - n_others;

% find rows with missing data
 TF = ismissing(tbl);
 tbl(any(TF,2),:);
 
% remove rows with missing data
 tbl = rmmissing(tbl);
 
% take log on family income
 tbl.FamilyIncome = log(tbl.FamilyIncome);

 % LASSO
tbl2 = tbl;
tbl2 = removevars(tbl, {'FamilyIncome', 'Occupation', 'MaritalStatus', 'Race', 'BirthState'});
tbl2 = [array2table([dummyvar(categorical(tbl.Occupation)),...
                           dummyvar(categorical(tbl.MaritalStatus)),...
                           dummyvar(categorical(tbl.Race)),...
                           dummyvar(categorical(tbl.BirthState))],...
         'VariableNames', {'Occupation_Blue', 'Occupation_Mng', 'Occupation_Service', 'Occupation_White',...
                           'MaritalStatus_Divorced_or_Widowed','MaritalStatus_Married','MaritalStatus_Single',...
                           'Race_Asian','Race_Black','Race_Hispanic','Race_White',...
                           'BirthState_MidWest','BirthState_Foreign','BirthState_Northeast','BirthState_South','BirthState_West'}),...
         tbl2];
tbl2 = removevars(tbl2, {'Occupation_White', 'MaritalStatus_Single', 'Race_White', 'BirthState_MidWest'});

n_others = size(tbl2,2) - n_measures;
tbl3 = table2array(tbl2(:, n_others+1:end));
tbl3 = (tbl3 - mean(tbl3))./std(tbl3);
D = x2fx(tbl3, 'Interaction');
D(:,1) = [];
D = [table2array(tbl2(:, 1:n_others)), D];

[B, fitInfo] = lasso(D, tbl.FamilyIncome, 'CV', 10);


a = subplot(1,2,2);
lassoPlot(B, fitInfo, 'PlotType', 'CV', 'parent', a)
l = legend('show');
set(l, 'Interpreter', 'latex');
set(get(a,'title'), 'Interpreter', 'latex');
set(get(a,'title'), 'String', 'Cross-Validated MSE of Lasso Fit (Female)');
set(get(a,'xlabel'), 'Interpreter', 'latex');
set(get(a,'ylabel'), 'Interpreter', 'latex');
set(a,'TickLabelInterpreter', 'latex');
 
fprintf('==============================================\n');
fprintf('|  LASSO RESULT (FEMALE)                       |\n');
fprintf('----------------------------------------------\n');
fprintf('|    Lambda @ 1SE: %f\n', fitInfo.Lambda1SE);
fprintf('|    Lambda @ Min: %f\n', fitInfo.LambdaMinMSE);
fprintf('|    SE @ 1SE: %f\n', fitInfo.SE(fitInfo.Index1SE));
fprintf('|    SE @ Min: %f\n', fitInfo.SE(fitInfo.IndexMinMSE));
fprintf('|    MSE @ 1SE: %f\n', fitInfo.MSE(fitInfo.Index1SE));
fprintf('|    MSE @ Min: %f\n', fitInfo.MSE(fitInfo.IndexMinMSE));
fprintf('===============================================\n');

mdl_lasso = fitlm(D(:, find(B(:,fitInfo.IndexMinMSE))), tbl.FamilyIncome);


  %% Reporting Error in Height/Weight vs Income (male only) - Table 4
% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
ReportedBMI = caesar.ReportedWeight./(caesar.ReportedHeight*0.001).^2;
AgeSquared = caesar.Age.^2;
Exp = max(caesar.Age - caesar.Education - 6,0);
ExpSquared = Exp.^2;
CarAge = 2001-caesar.CarYear;
WeightError = caesar.ReportedWeight - caesar.Weight;
BMIError = ReportedBMI - BMI;
HeightError = caesar.ReportedHeight - caesar.Stature;

tbl = [caesar, array2table([BMI, AgeSquared, Exp, ExpSquared, CarAge, HeightError, WeightError, BMIError], 'VariableNames', {'BMI', 'AgeSquared','Experience', 'ExperienceSquared', 'CarAge', 'HeightError', 'WeightError', 'BMIError'})];

% creating a subset that satisfies criteria (ex: selecting males, white collar)
 rows = strcmp(tbl.Gender, 'Male');
 tbl = tbl(rows,:);
 
% selecting variables for analysis
vars = {'FamilyIncome', 'Age', 'AgeSquared',  'Occupation', 'Education', ...
    'Fitness', 'MaritalStatus', 'Race', 'BirthState', 'WeightError', 'Weight'};
 
 tbl = tbl(:, vars);
 
% find rows with missing data
 TF = ismissing(tbl);
 tbl(any(TF,2),:);
 
% remove rows with missing data
 tbl = rmmissing(tbl);
 
% take log on family income
  tbl.FamilyIncome = log(tbl.FamilyIncome);
 
% fit a model
  mdl = fitlm(tbl, 'ResponseVar', 'WeightError')
  
coeff_estimate = mdl.Coefficients.Estimate';

% Bootstrap
N = size(tbl.FamilyIncome,1);

MAX_ITER = 1000;

bootstrap_coeff = [];
for i=1:MAX_ITER

    idx = randi(height(tbl), N,1);

    bootstrap_tbl = tbl(idx,:);
    bootstrap_mdl = fitlm(bootstrap_tbl, 'ResponseVar', 'WeightError');
    if mdl.NumCoefficients == length(bootstrap_mdl.Coefficients.Estimate)
        bootstrap_coeff = [bootstrap_coeff, bootstrap_mdl.Coefficients.Estimate];
    end
end

bootstrap_estimate = mean(bootstrap_coeff');
bootstrap_SE = std(bootstrap_coeff');
t_stat = coeff_estimate./bootstrap_SE;
fprintf('Bootstrap\n');
fprintf('============================================\n');
fprintf('%20s\tEstimate\tStd\tt-Stat\n', 'Variable');
for i=1:length(bootstrap_estimate)
    str = mdl.CoefficientNames{i};
    if length(str) > 20
        str = [str(1:5), '...', str(end-10:end)];
    end
    fprintf('%20s\t%f\t%f\t%f\n', str, coeff_estimate(i), bootstrap_SE(i), t_stat(i));
end

  
  %% Reporting Error in Height/Weight vs Income (Female only) - Table 4
% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
ReportedBMI = caesar.ReportedWeight./(caesar.ReportedHeight*0.001).^2;
AgeSquared = caesar.Age.^2;
Exp = max(caesar.Age - caesar.Education - 6,0);
ExpSquared = Exp.^2;
CarAge = 2001-caesar.CarYear;
WeightError = caesar.ReportedWeight - caesar.Weight;
BMIError = ReportedBMI - BMI;
HeightError = caesar.ReportedHeight - caesar.Stature;

tbl = [caesar, array2table([BMI, AgeSquared, Exp, ExpSquared, CarAge, HeightError, WeightError, BMIError], 'VariableNames', {'BMI', 'AgeSquared','Experience', 'ExperienceSquared', 'CarAge', 'HeightError', 'WeightError', 'BMIError'})];

% creating a subset that satisfies criteria (ex: selecting males, white collar)
 rows = strcmp(tbl.Gender, 'Female');
 tbl = tbl(rows,:);
 
% selecting variables for analysis
 vars = {'FamilyIncome',  'Age', 'AgeSquared', 'Occupation', 'Education', ...
     'Fitness', 'MaritalStatus', 'Race', 'BirthState', 'WeightError', 'Weight'};
 
 tbl = tbl(:, vars);
 
% find rows with missing data
 TF = ismissing(tbl);
 tbl(any(TF,2),:);
 
% remove rows with missing data
 tbl = rmmissing(tbl);
 
% take log on family income
  tbl.FamilyIncome = log(tbl.FamilyIncome);
 
% fit a model
  mdl = fitlm(tbl, 'ResponseVar', 'WeightError')

coeff_estimate = mdl.Coefficients.Estimate';
  
% Bootstrap
N = size(tbl.FamilyIncome,1);

MAX_ITER = 1000;

bootstrap_coeff = [];
for i=1:MAX_ITER

    idx = randi(height(tbl), N,1);

    bootstrap_tbl = tbl(idx,:);
    bootstrap_mdl = fitlm(bootstrap_tbl, 'ResponseVar', 'WeightError');
    if mdl.NumCoefficients == length(bootstrap_mdl.Coefficients.Estimate)
        bootstrap_coeff = [bootstrap_coeff, bootstrap_mdl.Coefficients.Estimate];
    end
end

bootstrap_estimate = mean(bootstrap_coeff');
bootstrap_SE = std(bootstrap_coeff');
t_stat = coeff_estimate./bootstrap_SE;
fprintf('Bootstrap\n');
fprintf('============================================\n');
fprintf('%20s\tEstimate\tStd\tt-Stat\n', 'Variable');
for i=1:length(bootstrap_estimate)
    str = mdl.CoefficientNames{i};
    if length(str) > 20
        str = [str(1:5), '...', str(end-10:end)];
    end
    fprintf('%20s\t%f\t%f\t%f\n', str, coeff_estimate(i), bootstrap_SE(i), t_stat(i));
end



  %% Reported Height/Weight vs Income (male only) - Table 5
% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
AgeSquared = caesar.Age.^2;
Exp = max(caesar.Age - caesar.Education - 6,0);
ExpSquared = Exp.^2;
CarAge = 2001-caesar.CarYear;
HeightError = caesar.ReportedHeight - caesar.Stature;
WeightError = caesar.ReportedWeight - caesar.Weight;

tbl = [caesar, array2table([BMI, AgeSquared, Exp, ExpSquared, CarAge, HeightError, WeightError], 'VariableNames', {'BMI', 'AgeSquared','Experience', 'ExperienceSquared', 'CarAge', 'HeightError', 'WeightError'})];

% creating a subset that satisfies criteria (ex: selecting males, white collar)
 rows = strcmp(tbl.Gender, 'Male');
 tbl = tbl(rows,:);
 
% selecting variables for analysis
 vars = {'FamilyIncome',  'Experience', 'ExperienceSquared', 'Occupation', 'Education', ...
         'MaritalStatus', 'Race',  'NumberOfChildren', 'ReportedHeight', 'ReportedWeight'};
 
 tbl = tbl(:, vars);
 
% find rows with missing data
 TF = ismissing(tbl);
 tbl(any(TF,2),:);
 
% remove rows with missing data
 tbl = rmmissing(tbl);
 
% take log on family income
  tbl.FamilyIncome = log(tbl.FamilyIncome);
 
% fit a model
  mdl = fitlm(tbl, 'ResponseVar', 'FamilyIncome')
 
coeff_estimate = mdl.Coefficients.Estimate';

% Bootstrap
N = size(tbl.FamilyIncome,1);

MAX_ITER = 1000;

bootstrap_coeff = [];
for i=1:MAX_ITER

    idx = randi(height(tbl), N,1);

    bootstrap_tbl = tbl(idx,:);
    bootstrap_mdl = fitlm(bootstrap_tbl, 'ResponseVar', 'FamilyIncome');
    if mdl.NumCoefficients == length(bootstrap_mdl.Coefficients.Estimate)
        bootstrap_coeff = [bootstrap_coeff, bootstrap_mdl.Coefficients.Estimate];
    end
end

bootstrap_estimate = mean(bootstrap_coeff');
bootstrap_SE = std(bootstrap_coeff');
t_stat = coeff_estimate./bootstrap_SE;
fprintf('Bootstrap\n');
fprintf('============================================\n');
fprintf('%20s\tEstimate\tStd\tt-Stat\n', 'Variable');
for i=1:length(bootstrap_estimate)
    str = mdl.CoefficientNames{i};
    if length(str) > 20
        str = [str(1:5), '...', str(end-10:end)];
    end
    fprintf('%20s\t%f\t%f\t%f\n', str, coeff_estimate(i), bootstrap_SE(i), t_stat(i));
end



  %% Reported Height/Weight vs Income (Female only) - Table 5
% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
AgeSquared = caesar.Age.^2;
Exp = max(caesar.Age - caesar.Education - 6,0);
ExpSquared = Exp.^2;
CarAge = 2001-caesar.CarYear;
HeightError = caesar.ReportedHeight - caesar.Stature;
WeightError = caesar.ReportedWeight - caesar.Weight;

tbl = [caesar, array2table([BMI, AgeSquared, Exp, ExpSquared, CarAge, HeightError, WeightError], 'VariableNames', {'BMI', 'AgeSquared','Experience', 'ExperienceSquared', 'CarAge', 'HeightError', 'WeightError'})];

% creating a subset that satisfies criteria (ex: selecting males, white collar)
 rows = strcmp(tbl.Gender, 'Female');
 tbl = tbl(rows,:);
 
% selecting variables for analysis
 vars = {'FamilyIncome',  'Experience', 'ExperienceSquared', 'Occupation', 'Education', ...
         'MaritalStatus', 'Race',  'NumberOfChildren', 'ReportedHeight', 'ReportedWeight'};
 
 tbl = tbl(:, vars);
 
% find rows with missing data
 TF = ismissing(tbl);
 tbl(any(TF,2),:);
 
% remove rows with missing data
 tbl = rmmissing(tbl);
 
% take log on family income
  tbl.FamilyIncome = log(tbl.FamilyIncome);
 
% fit a model
  mdl = fitlm(tbl, 'ResponseVar', 'FamilyIncome')

coeff_estimate = mdl.Coefficients.Estimate';

% Bootstrap
N = size(tbl.FamilyIncome,1);

MAX_ITER = 1000;

bootstrap_coeff = [];
for i=1:MAX_ITER

    idx = randi(height(tbl), N,1);

    bootstrap_tbl = tbl(idx,:);
    bootstrap_mdl = fitlm(bootstrap_tbl, 'ResponseVar', 'FamilyIncome');
    if mdl.NumCoefficients == length(bootstrap_mdl.Coefficients.Estimate)
        bootstrap_coeff = [bootstrap_coeff, bootstrap_mdl.Coefficients.Estimate];
    end
end

bootstrap_estimate = mean(bootstrap_coeff');
bootstrap_SE = std(bootstrap_coeff');
t_stat = coeff_estimate./bootstrap_SE;
fprintf('Bootstrap\n');
fprintf('============================================\n');
fprintf('%20s\tEstimate\tStd\tt-Stat\n', 'Variable');
for i=1:length(bootstrap_estimate)
    str = mdl.CoefficientNames{i};
    if length(str) > 20
        str = [str(1:5), '...', str(end-10:end)];
    end
    fprintf('%20s\t%f\t%f\t%f\n', str, coeff_estimate(i), bootstrap_SE(i), t_stat(i));
end


 %% Height/Weight vs Income (male only) - Tables 6
% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
AgeSquared = caesar.Age.^2;
Exp = max(caesar.Age - caesar.Education - 6,0);
ExpSquared = Exp.^2;
CarAge = 2001-caesar.CarYear;

tbl = [caesar, array2table([BMI, AgeSquared, Exp, ExpSquared, CarAge], 'VariableNames', {'BMI', 'AgeSquared','Experience', 'ExperienceSquared', 'CarAge'})];

% creating a subset that satisfies criteria (ex: selecting males, white collar)
 rows = strcmp(tbl.Gender, 'Male');
 tbl = tbl(rows,:);
 
% selecting variables for analysis
  vars = {'FamilyIncome',  'Experience', 'ExperienceSquared', 'Occupation', 'Education', ...
         'MaritalStatus', 'Race',  'NumberOfChildren', 'Stature', 'Weight'};
 
     
 tbl = tbl(:, vars);
 
% find rows with missing data
 TF = ismissing(tbl);
 tbl(any(TF,2),:);
 
% remove rows with missing data
 tbl = rmmissing(tbl);
 
% take log on family income
 tbl.FamilyIncome = log(tbl.FamilyIncome);
 
% fit a model
 mdl = fitlm(tbl, 'ResponseVar', 'FamilyIncome')
 
coeff_estimate = mdl.Coefficients.Estimate';

% Bootstrap
N = size(tbl.FamilyIncome,1);

MAX_ITER = 1000;

bootstrap_coeff = [];
for i=1:MAX_ITER
    
    idx = randi(height(tbl), N,1);

    bootstrap_tbl = tbl(idx,:);
    bootstrap_mdl = fitlm(bootstrap_tbl, 'ResponseVar', 'FamilyIncome');
    if mdl.NumCoefficients == length(bootstrap_mdl.Coefficients.Estimate)
        bootstrap_coeff = [bootstrap_coeff, bootstrap_mdl.Coefficients.Estimate];
    end
end

bootstrap_estimate = mean(bootstrap_coeff');
bootstrap_SE = std(bootstrap_coeff');
t_stat = coeff_estimate./bootstrap_SE;
fprintf('Bootstrap\n');
fprintf('============================================\n');
fprintf('%20s\tEstimate\tStd\tt-Stat\n', 'Variable');
for i=1:length(bootstrap_estimate)
    str = mdl.CoefficientNames{i};
    if length(str) > 20
        str = [str(1:5), '...', str(end-10:end)];
    end
    fprintf('%20s\t%f\t%f\t%f\n', str, coeff_estimate(i), bootstrap_SE(i), t_stat(i));
end



%% Height/Weight vs Income (Female only) - Table 6
% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
AgeSquared = caesar.Age.^2;
Exp = max(caesar.Age - caesar.Education - 6,0);
ExpSquared = Exp.^2;
CarAge = 2001-caesar.CarYear;

tbl = [caesar, array2table([BMI, AgeSquared, Exp, ExpSquared, CarAge], 'VariableNames', {'BMI', 'AgeSquared','Experience', 'ExperienceSquared', 'CarAge'})];

% creating a subset that satisfies criteria (ex: selecting males, white collar)
 rows = strcmp(tbl.Gender, 'Female');
 tbl = tbl(rows,:);
 
% selecting variables for analysis
 vars = {'FamilyIncome',  'Experience', 'ExperienceSquared', 'Occupation', 'Education', ...
         'MaritalStatus', 'Race',  'NumberOfChildren', 'Stature', 'Weight'};
 
 tbl = tbl(:, vars);
 
% find rows with missing data
 TF = ismissing(tbl);
 tbl(any(TF,2),:);
 
% remove rows with missing data
 tbl = rmmissing(tbl);
 
% take log on family income
 tbl.FamilyIncome = log(tbl.FamilyIncome);
 
% fit a model
 mdl = fitlm(tbl, 'ResponseVar', 'FamilyIncome')

coeff_estimate = mdl.Coefficients.Estimate';

% Bootstrap
N = size(tbl.FamilyIncome,1);

MAX_ITER = 1000;

bootstrap_coeff = [];
for i=1:MAX_ITER

    idx = randi(height(tbl), N,1);

    bootstrap_tbl = tbl(idx,:);
    bootstrap_mdl = fitlm(bootstrap_tbl, 'ResponseVar', 'FamilyIncome');
    if mdl.NumCoefficients == length(bootstrap_mdl.Coefficients.Estimate)
        bootstrap_coeff = [bootstrap_coeff, bootstrap_mdl.Coefficients.Estimate];
    end
end

bootstrap_estimate = mean(bootstrap_coeff');
bootstrap_SE = std(bootstrap_coeff');
t_stat = coeff_estimate./bootstrap_SE;
fprintf('Bootstrap\n');
fprintf('============================================\n');
fprintf('%20s\tEstimate\tStd\tt-Stat\n', 'Variable');
for i=1:length(bootstrap_estimate)
    str = mdl.CoefficientNames{i};
    if length(str) > 20
        str = [str(1:5), '...', str(end-10:end)];
    end
    fprintf('%20s\t%f\t%f\t%f\n', str, coeff_estimate(i), bootstrap_SE(i), t_stat(i));
end



%% Reported BMI vs Income (male only) - Table 7
% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
ReportedBMI = caesar.ReportedWeight./(caesar.ReportedHeight*0.001).^2;
AgeSquared = caesar.Age.^2;
Exp = max(caesar.Age - caesar.Education - 6,0);
ExpSquared = Exp.^2;
CarAge = 2001-caesar.CarYear;
WeightError = caesar.ReportedWeight - caesar.Weight;
BMIError = ReportedBMI - BMI;
HeightError = caesar.ReportedHeight - caesar.Stature;

tbl = [caesar, array2table([BMI, ReportedBMI, WeightError, BMIError, HeightError, AgeSquared,  Exp, ExpSquared, CarAge], 'VariableNames', {'BMI', 'ReportedBMI', 'WeightError', 'BMIError', 'HeightError', 'AgeSquared', 'Experience', 'ExperienceSquared', 'CarAge'})];

% creating a subset that satisfies criteria (ex: selecting males, white collar)
 rows = strcmp(tbl.Gender, 'Male');
 tbl = tbl(rows,:);
 
%selecting variables for analysis
   vars = {'FamilyIncome', 'Experience', 'ExperienceSquared', 'Occupation', 'Education', ...
           'MaritalStatus', 'Race', 'NumberOfChildren', 'ReportedBMI', 'ReportedHeight'};
  
 tbl = tbl(:, vars);
 
% find rows with missing data
 TF = ismissing(tbl);
 tbl(any(TF,2),:);
 
% remove rows with missing data
 tbl = rmmissing(tbl);
 
% take log on family income
 tbl.FamilyIncome = log(tbl.FamilyIncome);
 
% fit a model
 mdl = fitlm(tbl, 'ResponseVar', 'FamilyIncome')
 
coeff_estimate = mdl.Coefficients.Estimate';

% Bootstrap
N = size(tbl.FamilyIncome,1);

MAX_ITER = 1000;

bootstrap_coeff = [];
for i=1:MAX_ITER

    idx = randi(height(tbl), N,1);

    bootstrap_tbl = tbl(idx,:);
    bootstrap_mdl = fitlm(bootstrap_tbl, 'ResponseVar', 'FamilyIncome');
    if mdl.NumCoefficients == length(bootstrap_mdl.Coefficients.Estimate)
        bootstrap_coeff = [bootstrap_coeff, bootstrap_mdl.Coefficients.Estimate];
    end
end

bootstrap_estimate = mean(bootstrap_coeff');
bootstrap_SE = std(bootstrap_coeff');
t_stat = coeff_estimate./bootstrap_SE;
fprintf('Bootstrap\n');
fprintf('============================================\n');
fprintf('%20s\tEstimate\tStd\tt-Stat\n', 'Variable');
for i=1:length(bootstrap_estimate)
    str = mdl.CoefficientNames{i};
    if length(str) > 20
        str = [str(1:5), '...', str(end-10:end)];
    end
    fprintf('%20s\t%f\t%f\t%f\n', str, coeff_estimate(i), bootstrap_SE(i), t_stat(i));
end


%% Reported BMI vs Income (Female only) - Table 7
% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
ReportedBMI = caesar.ReportedWeight./(caesar.ReportedHeight*0.001).^2;
AgeSquared = caesar.Age.^2;
Exp = max(caesar.Age - caesar.Education - 6,0);
ExpSquared = Exp.^2;
CarAge = 2001-caesar.CarYear;
WeightError = caesar.ReportedWeight - caesar.Weight;
BMIError = ReportedBMI - BMI;
HeightError = caesar.ReportedHeight - caesar.Stature;

tbl = [caesar, array2table([BMI, ReportedBMI, WeightError, BMIError, HeightError, AgeSquared,  Exp, ExpSquared, CarAge], 'VariableNames', {'BMI', 'ReportedBMI', 'WeightError', 'BMIError', 'HeightError', 'AgeSquared', 'Experience', 'ExperienceSquared', 'CarAge'})];

% creating a subset that satisfies criteria (ex: selecting males, white collar)
 rows = strcmp(tbl.Gender, 'Female');
 tbl = tbl(rows,:);
 
%selecting variables for analysis
 vars = {'FamilyIncome', 'Experience', 'ExperienceSquared', 'Occupation', 'Education', ...
           'MaritalStatus', 'Race', 'NumberOfChildren', 'ReportedBMI', 'ReportedHeight'};
  
 tbl = tbl(:, vars);
 
% find rows with missing data
 TF = ismissing(tbl);
 tbl(any(TF,2),:);
 
% remove rows with missing data
 tbl = rmmissing(tbl);
 
% take log on family income
 tbl.FamilyIncome = log(tbl.FamilyIncome);
 
% fit a model
 mdl = fitlm(tbl, 'ResponseVar', 'FamilyIncome')

coeff_estimate = mdl.Coefficients.Estimate';

% Bootstrap
N = size(tbl.FamilyIncome,1);

MAX_ITER = 1000;

bootstrap_coeff = [];
for i=1:MAX_ITER

    idx = randi(height(tbl), N,1);

    bootstrap_tbl = tbl(idx,:);
    bootstrap_mdl = fitlm(bootstrap_tbl, 'ResponseVar', 'FamilyIncome');
    if mdl.NumCoefficients == length(bootstrap_mdl.Coefficients.Estimate)
        bootstrap_coeff = [bootstrap_coeff, bootstrap_mdl.Coefficients.Estimate];
    end
end

bootstrap_estimate = mean(bootstrap_coeff');
bootstrap_SE = std(bootstrap_coeff');
t_stat = coeff_estimate./bootstrap_SE;
fprintf('Bootstrap\n');
fprintf('============================================\n');
fprintf('%20s\tEstimate\tStd\tt-Stat\n', 'Variable');
for i=1:length(bootstrap_estimate)
    str = mdl.CoefficientNames{i};
    if length(str) > 20
        str = [str(1:5), '...', str(end-10:end)];
    end
    fprintf('%20s\t%f\t%f\t%f\n', str, coeff_estimate(i), bootstrap_SE(i), t_stat(i));
end


%% BMI vs Income (male only) - Tables 8 & 11 & 13
% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
AgeSquared = caesar.Age.^2;
Exp = max(caesar.Age - caesar.Education - 6,0);
ExpSquared = Exp.^2;
CarAge = 2001-caesar.CarYear;
tbl = [caesar, array2table([BMI, AgeSquared, Exp, ExpSquared, CarAge], 'VariableNames', {'BMI', 'AgeSquared', 'Experience', 'ExperienceSquared', 'CarAge'})];

% creating a subset that satisfies criteria (ex: selecting males, white collar)
 rows = strcmp(tbl.Gender, 'Male');
 tbl = tbl(rows,:);
 
%selecting variables for analysis

    vars = {'FamilyIncome',  'Experience', 'ExperienceSquared', 'Occupation', 'Education',  'MaritalStatus',  'Race', ...
         'NumberOfChildren', 'Fitness', 'CarModel',  'BirthState',  'Stature', 'BMI'};
     
  
 tbl = tbl(:, vars);
 
% find rows with missing data
 TF = ismissing(tbl);
 tbl(any(TF,2),:);
 
% remove rows with missing data
 tbl = rmmissing(tbl);
 
% take log on family income
 tbl.FamilyIncome = log(tbl.FamilyIncome);
 
% fit a model
 mdl = fitlm(tbl, 'ResponseVar', 'FamilyIncome')
 
coeff_estimate = mdl.Coefficients.Estimate';

% Bootstrap
N = size(tbl.FamilyIncome,1);

MAX_ITER = 1000;

bootstrap_coeff = [];
for i=1:MAX_ITER

    idx = randi(height(tbl), N,1);

    bootstrap_tbl = tbl(idx,:);
    bootstrap_mdl = fitlm(bootstrap_tbl, 'ResponseVar', 'FamilyIncome');
    if mdl.NumCoefficients == length(bootstrap_mdl.Coefficients.Estimate)
        bootstrap_coeff = [bootstrap_coeff, bootstrap_mdl.Coefficients.Estimate];
    end
end

bootstrap_estimate = mean(bootstrap_coeff');
bootstrap_SE = std(bootstrap_coeff');
t_stat = coeff_estimate./bootstrap_SE;
fprintf('Bootstrap\n');
fprintf('============================================\n');
fprintf('%20s\tEstimate\tStd\tt-Stat\n', 'Variable');
for i=1:length(bootstrap_estimate)
    str = mdl.CoefficientNames{i};
    if length(str) > 20
        str = [str(1:5), '...', str(end-10:end)];
    end
    fprintf('%20s\t%f\t%f\t%f\n', str, coeff_estimate(i), bootstrap_SE(i), t_stat(i));
end


 %% BMI vs Income (Female only) - Tables 8 & 11 & 13
% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
AgeSquared = caesar.Age.^2;
Exp = max(caesar.Age - caesar.Education - 6,0);
ExpSquared = Exp.^2;
CarAge = 2001-caesar.CarYear;
Hip2Waist = caesar.HipCircumference_Maximum./caesar.WaistCircumference_Pref.*100;
tbl = [caesar, array2table([BMI, AgeSquared, Exp, ExpSquared, CarAge, Hip2Waist], 'VariableNames', {'BMI', 'AgeSquared', 'Experience', 'ExperienceSquared', 'CarAge', 'Hip2Waist'})];

% creating a subset that satisfies criteria (ex: selecting males, white collar)
 rows = strcmp(tbl.Gender, 'Female');
 tbl = tbl(rows,:);
 
%selecting variables for analysis
   vars = {'FamilyIncome', 'Experience', 'ExperienceSquared', 'Occupation', 'Education', ...
           'MaritalStatus', 'Race', 'NumberOfChildren', 'Fitness',  'CarModel', 'BirthState', 'BMI', 'Stature', 'Hip2Waist'};
 
 tbl = tbl(:, vars);
 
% find rows with missing data
 TF = ismissing(tbl);
 tbl(any(TF,2),:);
 
% remove rows with missing data
 tbl = rmmissing(tbl);
 
% take log on family income
 tbl.FamilyIncome = log(tbl.FamilyIncome);
 
% fit a model
 mdl = fitlm(tbl, 'ResponseVar', 'FamilyIncome')

coeff_estimate = mdl.Coefficients.Estimate';

% Bootstrap
N = size(tbl.FamilyIncome,1);

MAX_ITER = 1000;

bootstrap_coeff = [];
for i=1:MAX_ITER
    
    idx = randi(height(tbl), N,1);

    bootstrap_tbl = tbl(idx,:);
    bootstrap_mdl = fitlm(bootstrap_tbl, 'ResponseVar', 'FamilyIncome');
    if mdl.NumCoefficients == length(bootstrap_mdl.Coefficients.Estimate)
        bootstrap_coeff = [bootstrap_coeff, bootstrap_mdl.Coefficients.Estimate];
    end
end

bootstrap_estimate = mean(bootstrap_coeff');
bootstrap_SE = std(bootstrap_coeff');
t_stat = coeff_estimate./bootstrap_SE;
fprintf('Bootstrap\n');
fprintf('============================================\n');
fprintf('%20s\tEstimate\tStd\tt-Stat\n', 'Variable');
for i=1:length(bootstrap_estimate)
    str = mdl.CoefficientNames{i};
    if length(str) > 20
        str = [str(1:5), '...', str(end-10:end)];
    end
    fprintf('%20s\t%f\t%f\t%f\n', str, coeff_estimate(i), bootstrap_SE(i), t_stat(i));
end


%% Various body measures vs Income - Table 9
% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
AgeSquared = caesar.Age.^2;
Exp = max(caesar.Age - caesar.Education - 6,0);
ExpSquared = Exp.^2;
CarAge = 2001-caesar.CarYear;

tbl = [caesar, array2table([BMI, AgeSquared, Exp, ExpSquared, CarAge], 'VariableNames', {'BMI', 'AgeSquared','Experience', 'ExperienceSquared', 'CarAge'})];

% creating a subset that satisfies criteria (ex: selecting males, white collar)
 gender = 'Male';
 rows = strcmp(tbl.Gender, gender);
 tbl = tbl(rows,:);
 
% selecting variables for analysis
vars = {'FamilyIncome',  'Experience', 'ExperienceSquared', 'Occupation', 'Education', ...
         'MaritalStatus', 'Race', 'NumberOfChildren'};
 

if strcmp(gender, 'Male') == 1
    vars = cat(2, vars, caesar.Properties.VariableNames(setdiff(26:70,[36,28,29,61,62,63])));
else
    vars = cat(2, vars, caesar.Properties.VariableNames(setdiff(26:70,[28,29,61,62,63])));
end
tbl = tbl(:, vars);

% find rows with missing data
 TF = ismissing(tbl);
 tbl(any(TF,2),:);
 
% remove rows with missing data
 tbl = rmmissing(tbl);
 
% take log on family income
 tbl.FamilyIncome = log(tbl.FamilyIncome);
 
% fit a model
 mdl = fitlm(tbl, 'ResponseVar', 'FamilyIncome') 

coeff_estimate = mdl.Coefficients.Estimate';

% Bootstrap
N = size(tbl.FamilyIncome,1);

MAX_ITER = 1000;

bootstrap_coeff = [];
for i=1:MAX_ITER

    idx = randi(height(tbl), N,1);

    bootstrap_tbl = tbl(idx,:);
    bootstrap_mdl = fitlm(bootstrap_tbl, 'ResponseVar', 'FamilyIncome');
    if mdl.NumCoefficients == length(bootstrap_mdl.Coefficients.Estimate)
        bootstrap_coeff = [bootstrap_coeff, bootstrap_mdl.Coefficients.Estimate];
    end
end

bootstrap_estimate = mean(bootstrap_coeff');
bootstrap_SE = std(bootstrap_coeff');
t_stat = coeff_estimate./bootstrap_SE;
fprintf('Bootstrap\n');
fprintf('============================================\n');
fprintf('%20s\tEstimate\tStd\tt-Stat\n', 'Variable');
for i=1:length(bootstrap_estimate)
    str = mdl.CoefficientNames{i};
    if length(str) > 20
        str = [str(1:5), '...', str(end-10:end)];
    end
    fprintf('%20s\t%f\t%f\t%f\n', str, coeff_estimate(i), bootstrap_SE(i), t_stat(i));
end

%% P3 vs Various body measures - Table 10
% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
Hip2Waist = caesar.HipCircumference_Maximum./caesar.WaistCircumference_Pref.*100;
tbl = [caesar, array2table([Hip2Waist], 'VariableNames', {'Hip2Waist'})];

% creating a subset that satisfies criteria (ex: selecting males, white collar)
gender = 'Female';
rows = strcmp(tbl.Gender, gender);
tbl = tbl(rows,:);
 
% selecting variables for analysis
vars = {'BraSize', 'CupSize', 'Var3'};
vars = cat(2, vars, caesar.Properties.VariableNames(setdiff(26:70,[28,29,61,62,63])));

tbl = tbl(:, vars);

% find rows with missing data
TF = ismissing(tbl);
tbl(any(TF,2),:);
 
% remove rows with missing data
tbl = rmmissing(tbl);
 
% fit a model
mdl = fitlm(tbl, 'ResponseVar', 'Var3') 

coeff_estimate = mdl.Coefficients.Estimate';

% Bootstrap
N = size(tbl.Var3,1);

MAX_ITER = 1000;

bootstrap_coeff = [];
for i=1:MAX_ITER
    %idx = randperm(height(tbl));
    %idx(N+1:end)=[];
    idx = randi(height(tbl), N,1);

    bootstrap_tbl = tbl(idx,:);
    bootstrap_mdl = fitlm(bootstrap_tbl, 'ResponseVar', 'Var3');
    if mdl.NumCoefficients == length(bootstrap_mdl.Coefficients.Estimate)
        bootstrap_coeff = [bootstrap_coeff, bootstrap_mdl.Coefficients.Estimate];
    end
end

bootstrap_estimate = mean(bootstrap_coeff');
bootstrap_SE = std(bootstrap_coeff');
t_stat = coeff_estimate./bootstrap_SE;
fprintf('Bootstrap\n');
fprintf('============================================\n');
fprintf('%20s\tEstimate\tStd\tt-Stat\n', 'Variable');
for i=1:length(bootstrap_estimate)
    str = mdl.CoefficientNames{i};
    if length(str) > 20
        str = [str(1:5), '...', str(end-10:end)];
    end
    fprintf('%20s\t%f\t%f\t%f\n', str, coeff_estimate(i), bootstrap_SE(i), t_stat(i));
end





 %% Body types from deep learning vs Income (male only) - Tables 12 & 14
% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
AgeSquared = caesar.Age.^2;
Exp = max(caesar.Age - caesar.Education - 6,0);
ExpSquared = Exp.^2;
CarAge = 2001-caesar.CarYear;
Interaction = caesar.Var1.*caesar.Var2;

tbl = [caesar, array2table([BMI, AgeSquared, Exp, ExpSquared, CarAge, Interaction], 'VariableNames', {'BMI', 'AgeSquared','Experience', 'ExperienceSquared', 'CarAge',  'Interaction'})];

% creating a subset that satisfies criteria (ex: selecting males, white collar)
  rows = strcmp(tbl.Gender, 'Male');
  tbl = tbl(rows,:);
 
%selecting variables for analysis
   vars = {'FamilyIncome',  'Experience', 'ExperienceSquared', 'Occupation', 'Education',  'MaritalStatus',  'Race', ...
         'NumberOfChildren', 'Fitness', 'CarModel',  'BirthState',  'Var1', 'Var2'};
    
%      vars = {'FamilyIncome',  'Experience', 'ExperienceSquared', 'Occupation', 'Education', 'MaritalStatus', 'Race',  ...
%          'NumberOfChildren', 'Var1', 'Var2'};


 tbl = tbl(:, vars);
 
% find rows with missing data
 TF = ismissing(tbl);
 tbl(any(TF,2),:);
 
% remove rows with missing data
 tbl = rmmissing(tbl);
 
% take log on family income
 tbl.FamilyIncome = log(tbl.FamilyIncome);

% fit a model
 mdl = fitlm(tbl, 'ResponseVar', 'FamilyIncome')

 
coeff_estimate = mdl.Coefficients.Estimate';

% Bootstrap
N = size(tbl.FamilyIncome,1);

MAX_ITER = 1000;

bootstrap_coeff = [];
for i=1:MAX_ITER

    idx = randi(height(tbl), N,1);

    bootstrap_tbl = tbl(idx,:);
    bootstrap_mdl = fitlm(bootstrap_tbl, 'ResponseVar', 'FamilyIncome');
    if mdl.NumCoefficients == length(bootstrap_mdl.Coefficients.Estimate)
        bootstrap_coeff = [bootstrap_coeff, bootstrap_mdl.Coefficients.Estimate];
    end
end

bootstrap_estimate = mean(bootstrap_coeff');
bootstrap_SE = std(bootstrap_coeff');
t_stat = coeff_estimate./bootstrap_SE;
fprintf('Bootstrap\n');
fprintf('============================================\n');
fprintf('%20s\tEstimate\tStd\tt-Stat\n', 'Variable');
for i=1:length(bootstrap_estimate)
    str = mdl.CoefficientNames{i};
    if length(str) > 20
        str = [str(1:5), '...', str(end-10:end)];
    end
    fprintf('%20s\t%f\t%f\t%f\n', str, coeff_estimate(i), bootstrap_SE(i), t_stat(i));
end



  %% Body types from deep learning vs Income (Female only) - Tables 12 & 14
% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
AgeSquared = caesar.Age.^2;
Exp = max(caesar.Age - caesar.Education - 6,0);
ExpSquared = Exp.^2;
CarAge = 2001-caesar.CarYear;
Interaction = caesar.Var1.*caesar.Var2;
SingleP3 = strcmp(caesar.MaritalStatus, 'Single').*caesar.Var3;


tbl = [caesar, array2table([BMI, AgeSquared, Exp, ExpSquared, CarAge, Interaction, SingleP3], 'VariableNames', {'BMI', 'AgeSquared','Experience', 'ExperienceSquared', 'CarAge',  'Interaction', 'SingleP3'})];

% creating a subset that satisfies criteria (ex: selecting males, white collar)

rows = strcmp(tbl.Gender, 'Female') ; 
tbl = tbl(rows,:);
 
%selecting variables for analysis
   vars = {'FamilyIncome',  'Experience', 'ExperienceSquared', 'Occupation', 'Education',  'MaritalStatus',  'Race', ...
        'NumberOfChildren', 'Fitness', 'CarModel',  'BirthState',   'Var1', 'Var2', 'Var3'};

%    vars = {'FamilyIncome',  'Experience', 'ExperienceSquared', 'Occupation', 'Education',  'MaritalStatus',  'Race', ...
%         'NumberOfChildren', 'Var1', 'Var2', 'Var3'};
    
        
 tbl = tbl(:, vars);
 
% find rows with missing data
 TF = ismissing(tbl);
 tbl(any(TF,2),:);
 
% remove rows with missing data
 tbl = rmmissing(tbl);
 
% take log on family income
 tbl.FamilyIncome = log(tbl.FamilyIncome);
 
% fit a model
 mdl = fitlm(tbl, 'ResponseVar', 'FamilyIncome')
 
 
coeff_estimate = mdl.Coefficients.Estimate';

% Bootstrap
N = size(tbl.FamilyIncome,1);
MAX_ITER = 1000;

bootstrap_coeff = [];
for i=1:MAX_ITER

    idx = randi(height(tbl), N,1);

    bootstrap_tbl = tbl(idx,:);
    bootstrap_mdl = fitlm(bootstrap_tbl, 'ResponseVar', 'FamilyIncome');
    if mdl.NumCoefficients == length(bootstrap_mdl.Coefficients.Estimate)
        bootstrap_coeff = [bootstrap_coeff, bootstrap_mdl.Coefficients.Estimate];
    end
end

bootstrap_estimate = mean(bootstrap_coeff');
bootstrap_SE = std(bootstrap_coeff');
t_stat = coeff_estimate./bootstrap_SE;
fprintf('Bootstrap\n');
fprintf('============================================\n');
fprintf('%20s\tEstimate\tStd\tt-Stat\n', 'Variable');
for i=1:length(bootstrap_estimate)
    str = mdl.CoefficientNames{i};
    if length(str) > 20
        str = [str(1:5), '...', str(end-10:end)];
    end
    fprintf('%20s\t%f\t%f\t%f\n', str, coeff_estimate(i), bootstrap_SE(i), t_stat(i));
end

%% IV: Body types from deep learning vs Income (Male only) - Tables 15 & 16

% Shoe Size Lookup
us2mm = containers.Map(...
{4,...
4.5,...
5,...
5.5,...
6,...
6.5,...
7,...
7.5,...
8,...
8.5,...
9,...
9.5,...
10,...
10.5,...
11,...
11.5,...
12,...
12.5,...
13,...
13.5,...
14},...
{220,...
224,...
229,...
233,...
237,...
241,...
245,...
250,...
254,...
258,...
262,...
267,...
271,...
275,...
279,...
283,...
288,...
292,...
296,...
300,...
305}...
);

% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
AgeSquared = caesar.Age.^2;
Exp = max(caesar.Age - caesar.Education - 6,0);
ExpSquared = Exp.^2;
CarAge = 2001-caesar.CarYear;
HeightError = caesar.ReportedHeight - caesar.Stature;

WeightError = caesar.ReportedWeight - caesar.Weight;

ShoeResidual = [];
for i=1:height(caesar)
    if isnan(caesar.ShoeSize(i))
        ShoeResidual(i,1) = NaN;
    else
        ShoeResidual(i,1) = us2mm(caesar.ShoeSize(i)) - caesar.FootLength(i);
    end
end
PantsResidual = caesar.PantsSizeWaist*25.4 - caesar.WaistCircumference_Pref;
JacketResidual = caesar.JacketSize*25.4 - caesar.ChestCircumference;

ShoeSize2 = caesar.ShoeSize.^2;
FootLength2 = caesar.FootLength.^2;

tbl = [caesar, array2table([BMI, AgeSquared, Exp, ExpSquared, CarAge, HeightError, WeightError,ShoeResidual,JacketResidual,PantsResidual, ShoeSize2, FootLength2], 'VariableNames', {'BMI', 'AgeSquared','Experience', 'ExperienceSquared', 'CarAge', 'HeightError', 'WeightError','ShoeResidual','JacketResidual','PantsResidual', 'ShoeSize2', 'FootLength2'})];


% creating a subset that satisfies criteria (ex: selecting males, white collar)
  rows = strcmp(tbl.Gender, 'Male');
  tbl = tbl(rows,:);
  
shoeres_model = fitlm(tbl.FootLength, tbl.ShoeSize);
tbl.ShoeResidual = shoeres_model.Residuals.Raw;

jackres_model = fitlm(tbl.ChestCircumference, tbl.JacketSize);
tbl.JacketResidual = jackres_model.Residuals.Raw;

pantres_model = fitlm(tbl.WaistCircumference_Pref, tbl.PantsSizeWaist);
tbl.PantResidual = pantres_model.Residuals.Raw;

%selecting variables for analysis
%use Stature and BMI in place of Var1 and Var2 for Table 15
    
     vars = {'FamilyIncome',  'Experience', 'ExperienceSquared', 'Occupation', 'Education',  'MaritalStatus',  'Race', ...
          'NumberOfChildren', 'Fitness', 'CarModel', 'BirthState', 'Site', 'ShoeResidual', 'JacketResidual', 'PantResidual', 'Var1', 'Var2'};


 tbl = tbl(:, vars);
 
% find rows with missing data
 TF = ismissing(tbl);
 tbl(any(TF,2),:);
 
% remove rows with missing data
 tbl = rmmissing(tbl);

 
vars = {'Experience', 'ExperienceSquared', 'Occupation', 'Education',  'MaritalStatus',  'Race', ...
         'NumberOfChildren', 'Fitness', 'CarModel', 'BirthState','Site',  'ShoeResidual', 'JacketResidual', 'PantResidual', 'Var1'};
tbl_v1 = tbl(:, vars);

% fit a model: 1st stage for Var1
 mdl_v1 = fitlm(tbl_v1, 'ResponseVar', 'Var1')

 coeff_estimate_v1 = mdl_v1.Coefficients.Estimate';
 
 % Bootstrap
N = size(tbl_v1.Var1,1);

MAX_ITER = 1000;

bootstrap_coeff_v1 = [];
for i=1:MAX_ITER

    idx = randi(height(tbl_v1), N,1);

    bootstrap_tbl_v1 = tbl_v1(idx,:);
    bootstrap_mdl_v1 = fitlm(bootstrap_tbl_v1, 'ResponseVar', 'Var1');
    if mdl_v1.NumCoefficients == length(bootstrap_mdl_v1.Coefficients.Estimate)
        bootstrap_coeff_v1 = [bootstrap_coeff_v1, bootstrap_mdl_v1.Coefficients.Estimate];
    end
end

bootstrap_estimate_v1 = mean(bootstrap_coeff_v1');
bootstrap_SE_v1 = std(bootstrap_coeff_v1');
t_stat_v1 = coeff_estimate_v1./bootstrap_SE_v1;
fprintf('Bootstrap\n');
fprintf('============================================\n');
fprintf('%20s\tEstimate\tStd\tt-Stat\n', 'Variable');
for i=1:length(bootstrap_estimate_v1)
    str = mdl_v1.CoefficientNames{i};
    if length(str) > 20
        str = [str(1:5), '...', str(end-10:end)];
    end
    fprintf('%20s\t%f\t%f\t%f\n', str, coeff_estimate_v1(i), bootstrap_SE_v1(i), t_stat_v1(i));
end


% define control function for Var1
CF1 = mdl_v1.Residuals(:,1);


% take log on family income
 tbl.FamilyIncome = log(tbl.FamilyIncome);
 
 vars = {'FamilyIncome',  'Experience', 'ExperienceSquared', 'Occupation', 'Education', 'MaritalStatus', 'Race',  ...
          'NumberOfChildren', 'Fitness', 'CarModel', 'BirthState',  'Site', 'Var1',   'Var2'};

     tbl2 = tbl(:, vars);
     
 tbl2 = [tbl2, CF1];
 tbl2.Properties.VariableNames{'Raw'} = 'ControlFunction1';
 
 % fit a model: 2nd stage
 mdl2 = fitlm(tbl2, 'ResponseVar', 'FamilyIncome')
 
coeff_estimate2 = mdl2.Coefficients.Estimate';

% Bootstrap
N = size(tbl2.FamilyIncome,1);

MAX_ITER = 1000;

bootstrap_coeff2 = [];
for i=1:MAX_ITER

    idx = randi(height(tbl2), N,1);

    bootstrap_tbl2 = tbl2(idx,:);
    bootstrap_mdl2 = fitlm(bootstrap_tbl2, 'ResponseVar', 'FamilyIncome');
    if mdl2.NumCoefficients == length(bootstrap_mdl2.Coefficients.Estimate)
        bootstrap_coeff2 = [bootstrap_coeff2, bootstrap_mdl2.Coefficients.Estimate];
    end
end

bootstrap_estimate2 = mean(bootstrap_coeff2');
bootstrap_SE2 = std(bootstrap_coeff2');
t_stat2 = coeff_estimate2./bootstrap_SE2;
fprintf('Bootstrap\n');
fprintf('============================================\n');
fprintf('%20s\tEstimate\tStd\tt-Stat\n', 'Variable');
for i=1:length(bootstrap_estimate2)
    str = mdl2.CoefficientNames{i};
    if length(str) > 20
        str = [str(1:5), '...', str(end-10:end)];
    end
    fprintf('%20s\t%f\t%f\t%f\n', str, coeff_estimate2(i), bootstrap_SE2(i), t_stat2(i));
end



 %% IV: Body types from deep learning vs Income (Female only) - Tables 15 & 16
 
% Shoe Size Lookup
us2mm = containers.Map(...
{4,...
4.5,...
5,...
5.5,...
6,...
6.5,...
7,...
7.5,...
8,...
8.5,...
9,...
9.5,...
10,...
10.5,...
11,...
11.5,...
12,...
12.5,...
13,...
13.5,...
14},...
{212,...
216,...
220,...
224,...
229,...
233,...
237,...
241,...
245,...
250,...
254,...
258,...
262,...
267,...
271,...
275,...
279,...
0,...
0,...
0,...
0}...
);

% adding new variables
BMI = caesar.Weight./(caesar.Stature*0.001).^2;
AgeSquared = caesar.Age.^2;
Exp = max(caesar.Age - caesar.Education - 6,0);
ExpSquared = Exp.^2;
CarAge = 2001-caesar.CarYear;
HeightError = caesar.ReportedHeight - caesar.Stature;
WeightError = caesar.ReportedWeight - caesar.Weight;
Hip2Waist = caesar.HipCircumference_Maximum./caesar.WaistCircumference_Pref.*100;

ShoeResidual = [];
for i=1:height(caesar)
    if isnan(caesar.ShoeSize(i))
        ShoeResidual(i,1) = NaN;
    else
        ShoeResidual(i,1) = us2mm(caesar.ShoeSize(i)) - caesar.FootLength(i);
    end
end
PantsResidual = caesar.PantsSizeWoman*25.4 - caesar.WaistCircumference_Pref;
BlouseResidual = caesar.BlouseSize*25.4 - caesar.ChestCircumference;
ShoeSize2 = caesar.ShoeSize.^2;
FootLength2 = caesar.FootLength.^2;

tbl = [caesar, array2table([BMI, AgeSquared, Exp, ExpSquared, CarAge, HeightError, WeightError,ShoeResidual,BlouseResidual,PantsResidual, ShoeSize2, FootLength2, Hip2Waist], 'VariableNames', {'BMI', 'AgeSquared','Experience', 'ExperienceSquared', 'CarAge', 'HeightError', 'WeightError','ShoeResidual','BlouseResidual','PantsResidual', 'ShoeSize2', 'FootLength2', 'Hip2Waist'})];


% creating a subset that satisfies criteria (ex: selecting males, white collar)
  rows = strcmp(tbl.Gender, 'Female');
  tbl = tbl(rows,:);
 
shoeres_model = fitlm(tbl.FootLength, tbl.ShoeSize);
tbl.ShoeResidual = shoeres_model.Residuals.Raw;

bloures_model = fitlm(tbl.ChestCircumference, tbl.BlouseSize);
tbl.BlouResidual = bloures_model.Residuals.Raw;

pantres_model = fitlm(tbl.WaistCircumference_Pref, tbl.PantsSizeWoman);
tbl.PantResidual = pantres_model.Residuals.Raw;

%selecting variables for analysis
%use Stature, BMI, and Hip2Waist in place of Var1, Var2, and Var3 for Table 15
     vars = {'FamilyIncome',  'Experience', 'ExperienceSquared', 'Occupation', 'Education', 'MaritalStatus', 'Race',  ...
             'NumberOfChildren', 'Fitness', 'CarModel', 'BirthState',  'ShoeResidual', 'BlouResidual', 'PantResidual', 'Var1', 'Var2', 'Var3'};

         
 tbl = tbl(:, vars);
 
% find rows with missing data
 TF = ismissing(tbl);
 tbl(any(TF,2),:);
 
% remove rows with missing data
 tbl = rmmissing(tbl);

vars = {'Experience', 'ExperienceSquared', 'Occupation', 'Education', 'MaritalStatus', 'Race',  ...
         'NumberOfChildren', 'Fitness', 'CarModel', 'BirthState',    'ShoeResidual', 'BlouResidual', 'PantResidual', 'Var1'};
tbl_v1 = tbl(:, vars);

% fit a model: 1st stage for Var1
 mdl_v1 = fitlm(tbl_v1, 'ResponseVar', 'Var1')

 coeff_estimate_v1 = mdl_v1.Coefficients.Estimate';
 
 % Bootstrap
N = size(tbl_v1.Var1,1);

MAX_ITER = 1000;

bootstrap_coeff_v1 = [];
for i=1:MAX_ITER

    idx = randi(height(tbl_v1), N,1);

    bootstrap_tbl_v1 = tbl_v1(idx,:);
    bootstrap_mdl_v1 = fitlm(bootstrap_tbl_v1, 'ResponseVar', 'Var1');
    if mdl_v1.NumCoefficients == length(bootstrap_mdl_v1.Coefficients.Estimate)
        bootstrap_coeff_v1 = [bootstrap_coeff_v1, bootstrap_mdl_v1.Coefficients.Estimate];
    end
end

bootstrap_estimate_v1 = mean(bootstrap_coeff_v1');
bootstrap_SE_v1 = std(bootstrap_coeff_v1');
t_stat_v1 = coeff_estimate_v1./bootstrap_SE_v1;
fprintf('Bootstrap\n');
fprintf('============================================\n');
fprintf('%20s\tEstimate\tStd\tt-Stat\n', 'Variable');
for i=1:length(bootstrap_estimate_v1)
    str = mdl_v1.CoefficientNames{i};
    if length(str) > 20
        str = [str(1:5), '...', str(end-10:end)];
    end
    fprintf('%20s\t%f\t%f\t%f\n', str, coeff_estimate_v1(i), bootstrap_SE_v1(i), t_stat_v1(i));
end


% define control function for Var1
CF1 = mdl_v1.Residuals(:,1);

 
% take log on family income
 tbl.FamilyIncome = log(tbl.FamilyIncome);
 
 vars = {'FamilyIncome',  'Experience', 'ExperienceSquared', 'Occupation', 'Education',  'MaritalStatus',  'Race', ...
         'NumberOfChildren', 'Fitness', 'CarModel', 'BirthState',  'Var1', 'Var2', 'Var3'};

     tbl2 = tbl(:, vars);
     
 tbl2 = [tbl2, CF1];
 tbl2.Properties.VariableNames{'Raw'} = 'ControlFunction1';
 
 % fit a model: 2nd stage
 mdl2 = fitlm(tbl2, 'ResponseVar', 'FamilyIncome')
 
coeff_estimate2 = mdl2.Coefficients.Estimate';

% Bootstrap
N = size(tbl2.FamilyIncome,1);

MAX_ITER = 1000;

bootstrap_coeff2 = [];
for i=1:MAX_ITER

    idx = randi(height(tbl2), N,1);

    bootstrap_tbl2 = tbl2(idx,:);
    bootstrap_mdl2 = fitlm(bootstrap_tbl2, 'ResponseVar', 'FamilyIncome');
    if mdl2.NumCoefficients == length(bootstrap_mdl2.Coefficients.Estimate)
        bootstrap_coeff2 = [bootstrap_coeff2, bootstrap_mdl2.Coefficients.Estimate];
    end
end

bootstrap_estimate2 = mean(bootstrap_coeff2');
bootstrap_SE2 = std(bootstrap_coeff2');
t_stat2 = coeff_estimate2./bootstrap_SE2;
fprintf('Bootstrap\n');
fprintf('============================================\n');
fprintf('%20s\tEstimate\tStd\tt-Stat\n', 'Variable');
for i=1:length(bootstrap_estimate2)
    str = mdl2.CoefficientNames{i};
    if length(str) > 20
        str = [str(1:5), '...', str(end-10:end)];
    end
    fprintf('%20s\t%f\t%f\t%f\n', str, coeff_estimate2(i), bootstrap_SE2(i), t_stat2(i));
end



