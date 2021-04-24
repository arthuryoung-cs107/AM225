fig1 = figure('Name', 'l2 error vs. N', 'Renderer', 'painters', 'Position', fig_pos(1, :));
xlabel('N_{fea}')
ylabel('l_2 error')
hold on

fig2 = figure('Name', 'l2 error vs. h', 'Renderer', 'painters', 'Position', fig_pos(2, :));
xlabel('h')
ylabel('l_2 error')
hold on


v_map = @(x, y) x*( 1 - 0.5* (y*y) ).^(0.5) ;
w_map = @(x, y) y*( 1 - 0.5* (x*x) ).^(0.5) ;
u_sol_xy_fun = @(x, y) (1 - v_map(x, y)^2 - w_map(x, y)^2)*exp(- v_map(x, y) );

l2_error = [];

N_test=500;

N_min = 10;
delta = 5;
N_max = 200;

for N=N_min:delta:N_max
  prob2_x = aysml_read(['../dat_dir/prob2_Ritz_Galerk_x_N', num2str(N), '_Ntest', num2str(N_test)]);
  prob2_y = aysml_read(['../dat_dir/prob2_Ritz_Galerk_y_N', num2str(N), '_Ntest', num2str(N_test)]);
  prob2_v = aysml_read(['../dat_dir/prob2_Ritz_Galerk_v_N', num2str(N), '_Ntest', num2str(N_test)]);
  prob2_w = aysml_read(['../dat_dir/prob2_Ritz_Galerk_w_N', num2str(N), '_Ntest', num2str(N_test)]);
  prob2_u = aysml_read(['../dat_dir/prob2_Ritz_Galerk_u_N', num2str(N), '_Ntest', num2str(N_test)]);

  u_sol = zeros(size(prob2_u));

  for i=1:size(u_sol, 1)
    for j=1:size(u_sol, 2)
      u_sol(i, j) = u_sol_xy_fun( prob2_x(i, j), prob2_y(i, j) );
    end
  end

  l2_error = [l2_error, (norm(u_sol-prob2_u ,'fro'))/(size(u_sol, 1)*size(u_sol, 2)) ];
end

figure(fig1.Number)
plot(N_min:delta:N_max, l2_error, ' o', 'Color', blue5, 'LineWidth', 1.5, 'DisplayName', '||u_{sol}-u_{fea}||/N_t^2')
legend('Show')

figure(fig2.Number)
plot( (N_min:delta:N_max).^(-1), l2_error, ' o', 'Color', red4, 'LineWidth', 1.5, 'DisplayName', '||u_{sol}-u_{fea}||/N_t^2')
legend('Show')
