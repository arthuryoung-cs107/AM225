fig1 = figure('Name', 'average l2 error vs. number of elements', 'Renderer', 'painters', 'Position', fig_pos(1, :));
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('N')
ylabel('l_2 error')
hold on

fig2 = figure('Name', 'average l2 error vs. h', 'Renderer', 'painters', 'Position', fig_pos(2, :));
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('h')
ylabel('l_2 error')
hold on

fig3 = figure('Name', 'average l2 error vs. number of elements', 'Renderer', 'painters', 'Position', fig_pos(3, :));
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('time elapsed (seconds)')
ylabel('l_2 error')
hold on

fig4 = figure('Name', 'x sol', 'Renderer', 'painters', 'Position', fig_pos(4, :));

N_test=1000;

N_min = 10;
delta = 5;
N_max = 1000;

u_sol_fun = @(x) exp(1-x).*sin(5*pi*x);

C0_data = aysml_read('../dat_dir/prob3_cube_C0_error');
C1_data = aysml_read('../dat_dir/prob3_altcube_C1_times');
C2_data = aysml_read('../dat_dir/prob3_altcube_C2_times');

C1_err = zeros(1, size(C1_data, 1));
C2_err = zeros(1, size(C2_data, 1));
% C1_xvals = zeros(1, size(C1_data, 1));
% C2_xvals = zeros(1, size(C2_data, 1));

i =  1;

for N=N_min:delta:N_max
  clf(fig4)
  C1_sol = aysml_read(['../dat_dir/prob3_altcube_C1_N', num2str(N),'_Ntest', num2str(N_test)]);
  C2_sol = aysml_read(['../dat_dir/prob3_altcube_C2_N', num2str(N),'_Ntest', num2str(N_test)]);

  C1_err(i) = norm(C1_sol(:, 2)-u_sol_fun(C1_sol(:, 1)) );
  C2_err(i) = norm(C2_sol(:, 2)-u_sol_fun(C2_sol(:, 1)) );

  C1_xvals = aysml_read(['../dat_dir/prob3_altcube_C1_N', num2str(N)]);
  C2_xvals = aysml_read(['../dat_dir/prob3_altcube_C2_N', num2str(N)]);

  figure(fig4.Number)
  xlabel('x')
  ylabel('u_{coeff}')
  hold on
  plot(linspace(1, 2, length(C1_xvals) ), C1_xvals, ' o', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', 'C1')
  plot(linspace(1, 2, length(C2_xvals) ), C2_xvals, ' o', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', 'C2')
  plot(linspace(1, 2, length(C1_xvals) ), u_sol_fun(linspace(1, 2, length(C1_xvals) )), ' o', 'Color', [0 0 0], 'LineWidth', 1.5, 'DisplayName', 'analytical')

  pause(1);

  i = i +1;
end

figure(fig1.Number)
plot(C1_data(:, 1), (C1_err')./C1_data(:, 1), ' o', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', 'C1 FEA')
plot(C2_data(:, 1), (C2_err')./C2_data(:, 1), ' o', 'Color', green4, 'LineWidth', 1.5, 'DisplayName', 'C2 FEA')
plot(C0_data(:, 1), C0_data(:, 4)./C0_data(:, 1), ' o', 'Color', [0 0 0], 'LineWidth', 1.5, 'DisplayName', 'C0 FEA')
legend('Show')

figure(fig2.Number)
plot(C1_data(:, 2), (C1_err')./C1_data(:, 1), ' o', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', 'C1 FEA')
plot(C2_data(:, 2), (C2_err')./C2_data(:, 1), ' o', 'Color', green4, 'LineWidth', 1.5, 'DisplayName', 'C2 FEA')
plot(C0_data(:, 2), C0_data(:, 4)./C0_data(:, 1), ' o', 'Color', [0 0 0], 'LineWidth', 1.5, 'DisplayName', 'C0 FEA')
legend('Show')

figure(fig3.Number)
plot(C1_data(:, 3), (C1_err')./C1_data(:, 1), ' o', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', 'C1 FEA')
plot(C2_data(:, 3), (C2_err')./C2_data(:, 1), ' o', 'Color', green4, 'LineWidth', 1.5, 'DisplayName', 'C2 FEA')
plot(C0_data(:, 3), C0_data(:, 4)./C0_data(:, 1), ' o', 'Color', [0 0 0], 'LineWidth', 1.5, 'DisplayName', 'C0 FEA')
legend('Show')
