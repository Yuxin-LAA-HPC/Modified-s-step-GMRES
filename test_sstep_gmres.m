%% Sparse example 1.
% load('/Users/yxma/Downloads/494_bus.mat')
% n = 494;
% A = Problem.A;

%% Sparse example 2.
load('/Users/yxma/Downloads/fs_183_6.mat')
n = 183;


%% Sparse example 3.
% load('/Users/yxma/Downloads/sherman2.mat')
% n = size(A, 1);

%% Modified s-step GMRES vs classical s-step GMRES.
b = ones(n, 1);
x0 = zeros(n, 1);
tolres = eps*n;
tolH = eps*sqrt(n);

normestA = norm(A, 'fro');
normb = norm(b);
smax = 32;

res = zeros(smax, 1);
res_tolH = zeros(smax, 1);
res_extraQR = zeros(smax, 1);
its = zeros(smax, 1);
its_tolH = zeros(smax, 1);
its_extraQR = zeros(smax, 1);
res_new = zeros(smax, 1);
res_tolH_new = zeros(smax, 1);
res_extraQR_new = zeros(smax, 1);
its_new = zeros(smax, 1);
its_tolH_new = zeros(smax, 1);
its_extraQR_new = zeros(smax, 1);
res_che = zeros(smax, 1);
res_tolH_che = zeros(smax, 1);
res_extraQR_che = zeros(smax, 1);
its_che = zeros(smax, 1);
its_tolH_che = zeros(smax, 1);
its_extraQR_che = zeros(smax, 1);
for s = 1:smax
    basis_info.type = 'newton';
    [x1, Z1, its1, flag1, error_res1, error_orth1] = gmres_sstep_norestart(A, x0, b, s, basis_info, tolres, tolH, 0, 0);
    res(s) = norm(A*x1 - b)/(normestA*norm(x1)+normb);
    its(s) = its1;
    [x2, Z2, its2, flag2, error_res2, error_orth2] = gmres_sstep_modified(A, x0, b, s, basis_info, tolres, tolH, 1);
    res_extraQR(s) = norm(A*x2 - b)/(normestA*norm(x2)+normb);
    its_extraQR(s) = its2;
    [x3, Z3, its3, flag3, error_res3, error_orth3] = gmres_sstep_norestart(A, x0, b, s, basis_info, tolres, tolH, 0, 1);
    res_tolH(s) = norm(A*x3 - b)/(normestA*norm(x3)+normb);
    its_tolH(s) = its3;

    basis_info.type = 'chebyshev';
    [x4, Z4, its4, flag4, error_res4, error_orth4] = gmres_sstep_norestart(A, x0, b, s, basis_info, tolres, tolH, 0, 0);
    res_che(s) = norm(A*x4 - b)/(normestA*norm(x4)+normb);
    its_che(s) = its4;
    [x5, Z5, its5, flag5, error_res5, error_orth5] = gmres_sstep_modified(A, x0, b, s, basis_info, tolres, tolH, 1);
    res_extraQR_che(s) = norm(A*x5 - b)/(normestA*norm(x5)+normb);
    its_extraQR_che(s) = its5;
    [x6, Z6, its6, flag6, error_res6, error_orth6] = gmres_sstep_norestart(A, x0, b, s, basis_info, tolres, tolH, 0, 1);
    res_tolH_che(s) = norm(A*x6 - b)/(normestA*norm(x6)+normb);
    its_tolH_che(s) = its6;
end

%% Draw figure to show that gmres_sstep_modified can employ much larger s.
figure
t = tiledlayout('flow','TileSpacing','compact');
nexttile
semilogy(res, 'LineWidth', 2, 'Marker','o');
hold on
semilogy(res_extraQR, 'LineWidth', 2, 'Marker','*');
hold on
semilogy(res_tolH, 'LineWidth', 2, 'Marker','+');
hold on
semilogy(res_che, 'LineWidth', 2, 'Marker','.');
hold on
semilogy(res_extraQR_che, 'LineWidth', 2, 'Marker','x');
hold on
semilogy(res_tolH_che, 'LineWidth', 2, 'Marker','square');
title('Relative backward error of final result');
xlabel('s');
ylabel('relative backward error');
set(gca,'FontSize', 18, 'FontWeight', 'normal')

nexttile
plot(its, 'LineWidth', 2, 'Marker', 'o');
hold on
plot(its_extraQR, 'LineWidth', 2, 'Marker', '*');
hold on
plot(its_tolH, 'LineWidth', 2, 'Marker', '+');
hold on
plot(its_che, 'LineWidth', 2, 'Marker', '.');
hold on
plot(its_extraQR_che, 'LineWidth', 2, 'Marker', 'x');
hold on
plot(its_tolH_che, 'LineWidth', 2, 'Marker', 'square');
title('Total number of iteratoins to converge');
xlabel('s');
ylabel('iter');
lgd = legend({['Classical s-step GMRES (Newton)'], ['Modified s-step GMRES with additional criteria (Newton)'], ['Classical s-step GMRES with additional criteria (Newton)'], ['Classical s-step GMRES (Chebyshev)'], ['Modified s-step GMRES with additional criteria (Chebyshev)'], ['Classical s-step GMRES with additional criteria (Chebyshev)']}, 'Location', 'southoutside');
lgd.Layout.Tile = 'south';
lgd.NumColumns = 2;
set(gca,'FontSize', 18, 'FontWeight', 'normal')
set(gcf,'Units','pixels','Position',[400 400 1000 450]);
