fig1 = figure('Name', 'cubic basis functions', 'Renderer', 'painters', 'Position', fig_pos(1, :));
xlabel('x')
ylabel('\phi (x) + i')
hold on

fig2 = figure('Name', 'cubic C1 basis functions', 'Renderer', 'painters', 'Position', fig_pos(2, :));
xlabel('x')
ylabel('\phi (x) + i')
hold on

fig3 = figure('Name', 'cubic C1 basis functions from C', 'Renderer', 'painters', 'Position', fig_pos(3, :));
xlabel('x')
ylabel('\phi (x) + i')
hold on

fig4 = figure('Name', 'cubic C1 basis functions from C', 'Renderer', 'painters', 'Position', fig_pos(4, :));
xlabel('x')
ylabel('d_x \phi (x) + i')
hold on

fig5 = figure('Name', 'cubic C2 basis ', 'Renderer', 'painters', 'Position', fig_pos(5, :));
xlabel('x')
ylabel(' \phi (x) + i')
hold on

fig6 = figure('Name', 'C2 phi check all ', 'Renderer', 'painters', 'Position', fig_pos(6, :));
xlabel('x')
ylabel(' \phi (x)')
hold on

fig7 = figure('Name', 'C2 grad phi check all ', 'Renderer', 'painters', 'Position', fig_pos(7, :));
xlabel('x')
ylabel('d \phi (x)')
hold on


omega = [1, 2];
N = 3;
N_full = 3*N + 1;
h = 1/(N_full);

x_vec = omega(1):(0.0001):(omega(2));

phi_mat = zeros(N_full, length(x_vec) );
phi_C1_mat = zeros(N+1, length(x_vec) );

phi0 = @(x) 1-0.75*(x.*x)+(0.25)*(x.*x.*x);

phiN= @(x) 1 + -x - (7.0/4.0)*(x.*x) - 0.5*(x.*x.*x);

figure(fig1.Number)
for i=1:N_full

  for j=1:length(x_vec)
      phi_mat(i, j) = phi_cubic(x_vec(j), i) + (i-1);
  end
  plot( x_vec, phi_mat(i, :), ' -', 'LineWidth', 1.5, 'DisplayName', ['\phi ', num2str(i -1)])
end
legend('Show')

figure(fig2.Number)
for i=1:(N+1)

  for j=1:length(x_vec)
      phi_C1_mat(i, j) = phi_cubic_C1(x_vec(j), i) + (i-1);
  end
  plot( x_vec, phi_C1_mat(i, :), ' -', 'LineWidth', 1.5, 'DisplayName', ['\phi ', num2str(i -1)])
end
legend('Show')


phi_C_check = aysml_read('../dat_dir/phi_check');
figure(fig3.Number)
for i=1:(size(phi_C_check, 2)-1)
  plot( phi_C_check(:, 1), phi_C_check(:, i+1), ' -', 'LineWidth', 1.5, 'DisplayName', ['\phi ', num2str(i -1)])
end
legend('Show')

grad_phi_C_check = aysml_read('../dat_dir/grad_phi_check');
figure(fig4.Number)
for i=1:(size(grad_phi_C_check, 2) -1 )
  plot( grad_phi_C_check(:, 1), grad_phi_C_check(:, i+1), ' -', 'LineWidth', 1.5, 'DisplayName', ['d_x \phi ', num2str(i -1)])
end
% legend('Show')

x_test = -2:0.01:2;
x_test0 = -1:0.01:2;
x_testN = -2:0.01:0;
phi_C2_vec = zeros(size(x_test));

for i=1:length(x_test)
  phi_C2_vec(i) = phi_cubic_C2(x_test(i));
end

figure(fig5.Number)
% plot( x_test, phi_C2_vec, ' -', 'LineWidth', 1.5, 'DisplayName', '')
% plot( x_test0, phi0(x_test0), ' -', 'LineWidth', 1.5, 'DisplayName', '')
plot( x_testN, phiN(x_testN), ' -', 'LineWidth', 1.5, 'DisplayName', '')


phi_C2_check_all = aysml_read('../dat_dir/phi_check_all');
figure(fig6.Number)
plot( phi_C2_check_all(:, 1), phi_C2_check_all(:, 2), ' -', 'LineWidth', 1.5, 'DisplayName', '')


figure(fig7.Number)
plot( phi_C2_check_all(:, 1), phi_C2_check_all(:, 3), ' -', 'LineWidth', 1.5, 'DisplayName', '')
