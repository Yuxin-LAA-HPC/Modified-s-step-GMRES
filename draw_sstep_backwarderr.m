%% Sparse example 1.
load('/Users/yxma/Downloads/494_bus.mat')
n = 494;
A = Problem.A;

%% Sparse example 2.
% load('/Users/yxma/Downloads/fs_183_6.mat')
% n = 183;

%% Sparse example 3.
% load('/Users/yxma/Downloads/sherman2.mat')
% n = size(A, 1);

%% sstep GMRES test.
b = ones(n, 1);
x0 = b;
tolres = eps*n;
tolH = eps*sqrt(n);

normestA = norm(A, 'fro');
normb = norm(b);
choose_s = [1, 6, 16];

figure
for i = 1:3
    s = choose_s(i);
    [x1, Z1, its1, flag1, error_res1] = gmres_sstep_norestart(A, x0, b, s, tolres, tolH, 1, 0);
    [x2, Z2, its2, flag2, error_res2] = gmres_sstep_modified(A, x0, b, s, tolres, tolH, 1);
    [x3, Z3, its3, flag3, error_res3] = gmres_sstep_norestart(A, x0, b, s, tolres, tolH, 1, 1);

    subplot(1, 3, i)
    semilogy(error_res1(1:its1), 'LineWidth', 2, 'Marker','o');
    hold on
    semilogy(error_res2(1:its2), 'LineWidth', 2, 'Marker','*');
    hold on
    semilogy(error_res3(1:its3), 'LineWidth', 2, 'Marker','+');
    hold on
    xlabel('iter');
    ylabel('backward error');
    set(gca,'FontSize', 18, 'FontWeight', 'normal')
end
legend({'Classical s-step GMRES', 'Modified s-step GMRES', ['Classical s-step GMRES' newline 'with additional criteria']}, 'Location', 'southeast');
set(gca,'FontSize', 18, 'FontWeight', 'normal')
set(gcf,'Units','pixels','Position',[400 400 1500 400]);