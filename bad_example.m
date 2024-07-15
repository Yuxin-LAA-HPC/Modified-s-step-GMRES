%% Bad dense example.
rng(1);
n = 20;
A = gallery('randsvd', [n, n], 1e5, 1);
[U, Sigma, V] = svd(A);
b = V(:, 4);
x0 = zeros(n, 1);
tolres = eps*n;
tolH = eps*sqrt(n);
s = 3;

[x1, Z1, its1, flag1, error_res1] = gmres_sstep_norestart(A, x0, b, s, tolres, tolH, 1, 0);
[x2, Z2, its2, flag2, error_res2] = gmres_sstep_modified(A, x0, b, s, tolres, tolH, 1);

figure
semilogy(error_res1(1:its1), 'LineWidth', 2, 'Marker','o');
hold on
semilogy(error_res2(1:its2), 'LineWidth', 2, 'Marker','*');
xlabel('iter');
ylabel('backward error');
set(gca,'FontSize', 18, 'FontWeight', 'normal')

