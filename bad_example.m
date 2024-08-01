%% Bad dense example.
rng(1);
n = 20;
A = gallery('randsvd', [n, n], 1e5, 1);
[U, Sigma, V] = svd(A);
b = V(:, 4);
x0 = zeros(n, 1);
tolres = eps*n;
tolH = eps*sqrt(n);

%% For s = 1
s = 1;
basis_info.type = 'monomial';
[x1, Z1, its1, flag1, error_res1, error_orth1, error_innerorth1] = gmres_sstep_norestart(A, x0, b, s, basis_info, tolres, tolH, 1, 0);
basis_info.type = 'newton';
[x2, Z2, its2, flag2, error_res2, error_orth2, error_innerorth2] = gmres_sstep_norestart(A, x0, b, s, basis_info, tolres, tolH, 1, 0);
basis_info.type = 'chebyshev';
[x3, Z3, its3, flag3, error_res3, error_orth3, error_innerorth3] = gmres_sstep_norestart(A, x0, b, s, basis_info, tolres, tolH, 1, 0);

figure
subplot(2, 3, 1);
semilogy(error_res1(1:its1), 'LineWidth', 2, 'Marker','o');
hold on
semilogy(error_res2(1:its2), 'LineWidth', 2, 'Marker','*');
hold on
semilogy(error_res3(1:its3), 'LineWidth', 2, 'Marker','+');
title('s = 1');
xlabel('iter');
ylabel('relative backward error');
ylim([1e-18, 1e-3]);
yticks([1e-16, 1e-12, 1e-8, 1e-4]);
set(gca,'FontSize', 18, 'FontWeight', 'normal')
legend('Monomial', 'Newton', 'Chebyshev', 'Location', 'southeast');

subplot(2, 3, 2);
semilogy(error_orth1(1:its1), 'LineWidth', 2, 'Marker','o');
hold on
semilogy(error_orth2(1:its2), 'LineWidth', 2, 'Marker','*');
hold on
semilogy(error_orth3(1:its3), 'LineWidth', 2, 'Marker','+');
title('s = 1');
xlabel('iter');
ylabel('\kappa(B_{1:iter}D_{1:iter})');
ylim([0.1, 1e15]);
yticks([1e0, 1e5, 1e10, 1e15])
set(gca,'FontSize', 18, 'FontWeight', 'normal')
legend('Monomial', 'Newton', 'Chebyshev');

subplot(2, 3, 3);
semilogy(error_innerorth1(1:ceil(its1/s)), 'LineWidth', 2, 'Marker','o');
hold on
semilogy(error_innerorth2(1:ceil(its2/s)), 'LineWidth', 2, 'Marker','*');
hold on
semilogy(error_innerorth3(1:ceil(its3/s)), 'LineWidth', 2, 'Marker','+');
title('s = 1');
xlabel('i');
ylabel('\kappa(B_{(i-1)s:is}D_{(i-1)s:is})');
ylim([0.1, 1e15]);
yticks([1e0, 1e5, 1e10, 1e15])
set(gca,'FontSize', 18, 'FontWeight', 'normal')
legend('Monomial', 'Newton', 'Chebyshev');

%% For s = 3
s = 3;
basis_info.type = 'monomial';
[x1, Z1, its1, flag1, error_res1, error_orth1, error_innerorth1] = gmres_sstep_norestart(A, x0, b, s, basis_info, tolres, tolH, 1, 0);
basis_info.type = 'newton';
[x2, Z2, its2, flag2, error_res2, error_orth2, error_innerorth2] = gmres_sstep_norestart(A, x0, b, s, basis_info, tolres, tolH, 1, 0);
basis_info.type = 'chebyshev';
[x3, Z3, its3, flag3, error_res3, error_orth3, error_innerorth3] = gmres_sstep_norestart(A, x0, b, s, basis_info, tolres, tolH, 1, 0);


subplot(2, 3, 4);
semilogy(error_res1(1:its1), 'LineWidth', 2, 'Marker','o');
hold on
semilogy(error_res2(1:its2), 'LineWidth', 2, 'Marker','*');
hold on
semilogy(error_res3(1:its3), 'LineWidth', 2, 'Marker','+');
title('s = 3');
xlabel('iter');
ylabel('relative backward error');
ylim([1e-18, 1e-3]);
yticks([1e-16, 1e-12, 1e-8, 1e-4]);
set(gca,'FontSize', 18, 'FontWeight', 'normal')
legend('Monomial', 'Newton', 'Chebyshev', 'Location', 'southeast');

subplot(2, 3, 5);
semilogy(error_orth1(1:its1), 'LineWidth', 2, 'Marker','o');
hold on
semilogy(error_orth2(1:its2), 'LineWidth', 2, 'Marker','*');
hold on
semilogy(error_orth3(1:its3), 'LineWidth', 2, 'Marker','+');
title('s = 3');
xlabel('iter');
ylabel('\kappa(B_{1:iter}D_{1:iter})');
ylim([0.1, 1e15]);
yticks([1e0, 1e5, 1e10, 1e15])
set(gca,'FontSize', 18, 'FontWeight', 'normal')
legend('Monomial', 'Newton', 'Chebyshev', 'Location', 'southeast');

subplot(2, 3, 6);
semilogy(error_innerorth1(1:ceil(its1/s)), 'LineWidth', 2, 'Marker','o');
hold on
semilogy(error_innerorth2(1:ceil(its2/s)), 'LineWidth', 2, 'Marker','*');
hold on
semilogy(error_innerorth3(1:ceil(its3/s)), 'LineWidth', 2, 'Marker','+');
title('s = 3');
xlabel('i');
ylabel('\kappa(B_{(i-1)s:is}D_{(i-1)s:is})');
ylim([0.1, 1e15]);
yticks([1e0, 1e5, 1e10, 1e15])
set(gca,'FontSize', 18, 'FontWeight', 'normal')
legend('Monomial', 'Newton', 'Chebyshev');
set(gcf,'Units','pixels','Position',[300 300 1500 600]);


