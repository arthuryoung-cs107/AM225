fig1 = figure('Name', 'C0 cubic basis functions', 'Renderer', 'painters', 'Position', fig_pos(1, :));
xlabel('x')
ylabel('\phi (x) + i')
hold on

fig2 = figure('Name', 'C1 cubic basis functions', 'Renderer', 'painters', 'Position', fig_pos(2, :));
xlabel('x')
ylabel('\phi (x) + i')
hold on

fig3 = figure('Name', 'C1 cubic basis gradient', 'Renderer', 'painters', 'Position', fig_pos(3, :));
xlabel('x')
ylabel('\phi (x) + i')
hold on

fig4 = figure('Name', 'C2 cubic basis functions', 'Renderer', 'painters', 'Position', fig_pos(4, :));
xlabel('x')
ylabel('\phi (x) + i')
hold on

fig5 = figure('Name', 'C2 cubic basis gradient', 'Renderer', 'painters', 'Position', fig_pos(5, :));
xlabel('x')
ylabel('\phi (x) + i')
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
legend('Show', 'Location', 'NorthWest')

phi_C1_check = aysml_read('../dat_dir/phi_check_C1');
figure(fig2.Number)
for i=1:(size(phi_C1_check, 2)-1)
  plot( phi_C1_check(:, 1), phi_C1_check(:, i+1), ' -', 'LineWidth', 1.5, 'DisplayName', ['\phi ', num2str(i -1)])
end
legend('Show', 'Location', 'NorthWest')

grad_phi_C1_check = aysml_read('../dat_dir/grad_phi_check_C1');
figure(fig3.Number)
for i=1:(size(grad_phi_C1_check, 2) -1 )
  plot( grad_phi_C1_check(:, 1), grad_phi_C1_check(:, i+1), ' -', 'LineWidth', 1.5, 'DisplayName', ['d_x \phi ', num2str(i -1)])
end
legend('Show', 'Location', 'SouthWest')

phi_C2_check = aysml_read('../dat_dir/phi_check_C2');
figure(fig4.Number)
for i=1:(size(phi_C2_check, 2)-1)
  plot( phi_C2_check(:, 1), phi_C2_check(:, i+1), ' -', 'LineWidth', 1.5, 'DisplayName', ['\phi ', num2str(i -1)])
end
legend('Show', 'Location', 'NorthWest')

grad_phi_C2_check = aysml_read('../dat_dir/grad_phi_check_C2');
figure(fig5.Number)
for i=1:(size(grad_phi_C2_check, 2) -1 )
  plot( grad_phi_C2_check(:, 1), grad_phi_C2_check(:, i+1), ' -', 'LineWidth', 1.5, 'DisplayName', ['d_x \phi ', num2str(i -1)])
end
legend('Show', 'Location', 'SouthEast')
