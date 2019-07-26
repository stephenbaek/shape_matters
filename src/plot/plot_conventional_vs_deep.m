function plot_conventional_vs_deep(tbl)
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
end