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

%% Modified s-step GMRES vs classical s-step GMRES.
b = ones(n, 1);
x0 = b;
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
for s = 1:smax
    [x1, Z1, its1, flag1, error_res1] = gmres_sstep_norestart(A, x0, b, s, tolres, tolH, 0, 0);
    res(s) = norm(A*x1 - b)/(normestA*norm(x1)+normb);
    its(s) = its1;
    [x2, Z2, its2, flag2, error_res2] = gmres_sstep_modified(A, x0, b, s, tolres, tolH, 0);
    res_extraQR(s) = norm(A*x2 - b)/(normestA*norm(x2)+normb);
    its_extraQR(s) = its2;
    [x3, Z3, its3, flag3, error_res3] = gmres_sstep_norestart(A, x0, b, s, tolres, tolH, 0, 1);
    res_tolH(s) = norm(A*x3 - b)/(normestA*norm(x3)+normb);
    its_tolH(s) = its3;
end

%% Draw figure to show that gmres_sstep_modified can employ much larger s.
figure
subplot(1, 2, 1)
semilogy(res, 'LineWidth', 2, 'Marker','o');
hold on
semilogy(res_extraQR, 'LineWidth', 2, 'Marker','*');
hold on
semilogy(res_tolH, 'LineWidth', 2, 'Marker','+');
xlabel('s');
ylabel('backward error');
set(gca,'FontSize', 18, 'FontWeight', 'normal')
subplot(1, 2, 2)
plot(its, 'LineWidth', 2, 'Marker', 'o');
hold on
plot(its_extraQR, 'LineWidth', 2, 'Marker', '*');
hold on
plot(its_tolH, 'LineWidth', 2, 'Marker', '+');
xlabel('s');
ylabel('iter');
legend({'Classical s-step GMRES', 'Modified s-step GMRES', ['Classical s-step GMRES' newline 'with additional criteria']}, 'Location', 'southeast');
set(gca,'FontSize', 18, 'FontWeight', 'normal')
set(gcf,'Units','pixels','Position',[400 400 1000 400]);