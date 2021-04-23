fig1 = figure('Name', 'u, Ritz Galerk, xy', 'Renderer', 'painters', 'Position', fig_pos(1, :));
xlabel('x')
ylabel('y')
zlabel('u')
view(40, 35)
hold on

fig2 = figure('Name', 'u, Ritz Galerk, vw', 'Renderer', 'painters', 'Position', fig_pos(2, :));
xlabel('v')
ylabel('w')
zlabel('u')
view(40, 35)
hold on

fig3 = figure('Name', 'u, analytical, xy', 'Renderer', 'painters', 'Position', fig_pos(3, :));
xlabel('x')
ylabel('y')
zlabel('u')
view(40, 35)
hold on

fig4 = figure('Name', 'u, analytical, wv', 'Renderer', 'painters', 'Position', fig_pos(4, :));
xlabel('w')
ylabel('v')
zlabel('u')
view(40, 35)
hold on

fig5 = figure('Name', 'l2 error', 'Renderer', 'painters', 'Position', fig_pos(5, :));
xlabel('N')
ylabel('l_2 error')
hold on


v_map = @(x, y) x*( 1 - 0.5* (y*y) ).^(0.5) ;
w_map = @(x, y) y*( 1 - 0.5* (x*x) ).^(0.5) ;
u_sol_xy_fun = @(x, y) (1 - v_map(x, y)^2 - w_map(x, y)^2)*exp(- v_map(x, y) );

l2_error = [];

N_max = 200;

for N=10:5:N_max
  prob2_x = aysml_read(['../dat_dir/prob2_Ritz_Galerk_x_N', num2str(N)]);
  prob2_y = aysml_read(['../dat_dir/prob2_Ritz_Galerk_y_N', num2str(N)]);
  prob2_v = aysml_read(['../dat_dir/prob2_Ritz_Galerk_v_N', num2str(N)]);
  prob2_w = aysml_read(['../dat_dir/prob2_Ritz_Galerk_w_N', num2str(N)]);
  prob2_u = aysml_read(['../dat_dir/prob2_Ritz_Galerk_u_N', num2str(N)]);

  u_sol = zeros(size(prob2_u));

  for i=1:size(u_sol, 1)
    for j=1:size(u_sol, 2)
      u_sol(i, j) = u_sol_xy_fun( prob2_x(i, j), prob2_y(i, j) );
    end
  end

  l2_error = [l2_error, (norm(u_sol-prob2_u ,'fro'))/(size(u_sol, 1)*size(u_sol, 2)) ];

  % figure(fig1.Number);
  % surf(prob2_x, prob2_y, prob2_u)
  %
  % figure(fig2.Number);
  % surf(prob2_v, prob2_w, prob2_u)

  % u_sol = zeros(size(prob2_u));
  %
  % for i=1:size(u_sol, 1)
  %   for j=1:size(u_sol, 2)
  %     u_sol(i, j) = u_sol_xy_fun( prob2_x(i, j), prob2_y(i, j) );
  %   end
  % end

  % figure(fig3.Number)
  % surf(prob2_x, prob2_y, u_sol)
  %
  % figure(fig4.Number)
  % surf(prob2_v, prob2_w, u_sol)

end

figure(fig5.Number)
plot(10:5:N_max, l2_error, ' o', 'Color', blue5, 'LineWidth', 1.5, 'DisplayName', '||u_sol-u_RG||/(N^2)')
