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

fig5 = figure('Name', 'source, xy', 'Renderer', 'painters', 'Position', fig_pos(5, :));
xlabel('x')
ylabel('y')
zlabel('f')
view(40, 35)
hold on

fig6 = figure('Name', 'source, wv', 'Renderer', 'painters', 'Position', fig_pos(6, :));
xlabel('w')
ylabel('v')
zlabel('f')
view(40, 35)
hold on


v_map = @(x, y) x*( 1 - 0.5* (y*y) ).^(0.5) ;
w_map = @(x, y) y*( 1 - 0.5* (x*x) ).^(0.5) ;
u_sol_xy_fun = @(x, y) (1 - v_map(x, y)^2 - w_map(x, y)^2)*exp(- v_map(x, y) );
u_source_fun = @(x, y) (3.0 + (v_map(x, y) - 4.0 )*v_map(x, y) + w_map(x, y)*w_map(x, y))*exp(- v_map(x, y) );

l2_error = [];

N_test=40;

N_min = 10;
delta = 5;
N_max = 100;

for N=N_min:delta:N_max
  clf(fig1)
  clf(fig2)
  clf(fig3)
  clf(fig4)
  clf(fig5)
  clf(fig6)

  prob2_x = aysml_read(['../dat_dir/prob2_Ritz_Galerk_x_N', num2str(N), '_Ntest', num2str(N_test)]);
  prob2_y = aysml_read(['../dat_dir/prob2_Ritz_Galerk_y_N', num2str(N), '_Ntest', num2str(N_test)]);
  prob2_v = aysml_read(['../dat_dir/prob2_Ritz_Galerk_v_N', num2str(N), '_Ntest', num2str(N_test)]);
  prob2_w = aysml_read(['../dat_dir/prob2_Ritz_Galerk_w_N', num2str(N), '_Ntest', num2str(N_test)]);
  prob2_u = aysml_read(['../dat_dir/prob2_Ritz_Galerk_u_N', num2str(N), '_Ntest', num2str(N_test)]);

  u_sol = zeros(size(prob2_u));
  u_source = zeros(size(prob2_u));

  for i=1:size(u_sol, 1)
    for j=1:size(u_sol, 2)
      u_sol(i, j) = u_sol_xy_fun( prob2_x(i, j), prob2_y(i, j) );
      u_source(i, j) = u_source_fun( prob2_x(i, j), prob2_y(i, j) );
    end
  end

  l2_error = [l2_error, (norm(u_sol-prob2_u ,'fro'))/(size(u_sol, 1)*size(u_sol, 2)) ];

  figure(fig1.Number);
  surf(prob2_x, prob2_y, prob2_u)

  figure(fig2.Number);
  surf(prob2_v, prob2_w, prob2_u)

  figure(fig3.Number)
  surf(prob2_x, prob2_y, u_sol)

  figure(fig4.Number)
  surf(prob2_v, prob2_w, u_sol)

  figure(fig5.Number)
  surf(prob2_x, prob2_y, u_source)

  figure(fig6.Number)
  surf(prob2_v, prob2_w, u_source)

  pause(1)

end
