%% Sparse example 1.
load('/Users/yxma/Downloads/494_bus.mat')
n = 494;
A = Problem.A;

%% Sparse example 2.
load('/Users/yxma/Downloads/fs_183_6.mat')
n = 183;

%% Sparse example 3.
% load('/Users/yxma/Downloads/sherman2.mat')
% n = size(A, 1);

%% sstep GMRES test.
b = ones(n, 1);
x0 = zeros(n, 1);%b;
tolres = eps*n;
tolH = eps*sqrt(n);

normestA = norm(A, 'fro');
normb = norm(b);
choose_s = [1, 4, 16];

figure
t = tiledlayout(1, 3);
for i = 1:size(choose_s, 2)
    s = choose_s(i);

    basis_info.type = 'newton';
    [x1, Z1, its1, flag1, error_res1, error_orth1] = gmres_sstep_norestart(A, x0, b, s, basis_info, tolres, tolH, 1, 0);
    [x2, Z2, its2, flag2, error_res2, error_orth2] = gmres_sstep_modified(A, x0, b, s, basis_info, tolres, tolH, 1);
    [x3, Z3, its3, flag3, error_res3, error_orth3] = gmres_sstep_norestart(A, x0, b, s, basis_info, tolres, tolH, 1, 1);

    basis_info.type = 'chebyshev';
    [x4, Z4, its4, flag4, error_res4, error_orth4] = gmres_sstep_norestart(A, x0, b, s, basis_info, tolres, tolH, 1, 0);
    [x5, Z5, its5, flag5, error_res5, error_orth5] = gmres_sstep_modified(A, x0, b, s, basis_info, tolres, tolH, 1);
    [x6, Z6, its6, flag6, error_res6, error_orth6] = gmres_sstep_norestart(A, x0, b, s, basis_info, tolres, tolH, 1, 1);

    nexttile
    semilogy(error_orth1(1:its1), 'LineWidth', 2, 'Marker','o');
    hold on
    semilogy(error_orth2(1:its2), 'LineWidth', 2, 'Marker','*');
    hold on
    semilogy(error_orth3(1:its3), 'LineWidth', 2, 'Marker','+');
    hold on
    semilogy(error_orth4(1:its4), 'LineWidth', 2, 'Marker','.');
    hold on
    semilogy(error_orth5(1:its5), 'LineWidth', 2, 'Marker','x');
    hold on
    semilogy(error_orth6(1:its6), 'LineWidth', 2, 'Marker','square');
    hold on
    xlabel('iter');
    ylabel('condition number of B');
    set(gca,'FontSize', 18, 'FontWeight', 'normal')
end
lgd = legend({'Classical s-step GMRES (Newton)', 'Modified s-step GMRES with additional criteria (Newton)', 'Classical s-step GMRES with additional criteria (Newton)', 'Classical s-step GMRES (Chebyshev)', 'Modified s-step GMRES with additional criteria (Chebyshev)', 'Classical s-step GMRES with additional criteria (Chebyshev)'}, 'Location', 'bestoutside');
set(gca,'FontSize', 18, 'FontWeight', 'normal')
lgd.NumColumns = 2;
lgd.Layout.Tile = 'south';
set(gca,'FontSize', 18, 'FontWeight', 'normal')
set(gcf,'Units','pixels','Position',[400 400 1200 400]);

