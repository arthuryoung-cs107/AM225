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
ylabel('grad \phi (x) + i')
hold on

omega = [1, 2];
N = 3;
N_full = 3*N + 1;
h = 1/(N_full);

x_vec = omega(1):(0.0001):(omega(2));

phi_mat = zeros(N_full, length(x_vec) );
phi_C1_mat = zeros(N+1, length(x_vec) );

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
for i=1:(size(phi_C_check, 2) -1 )
  plot( phi_C_check(:, 1), phi_C_check(:, i+1), ' -', 'LineWidth', 1.5, 'DisplayName', ['\phi ', num2str(i -1)])
end
legend('Show')

grad_phi_C_check = aysml_read('../dat_dir/grad_phi_check');
figure(fig4.Number)
for i=1:(size(grad_phi_C_check, 2) -1 )
  plot( grad_phi_C_check(:, 1), grad_phi_C_check(:, i+1), ' -', 'LineWidth', 1.5, 'DisplayName', ['d \phi ', num2str(i -1)])
end
legend('Show')
