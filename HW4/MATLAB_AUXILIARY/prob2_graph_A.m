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

prob2_x_20 = aysml_read('../dat_dir/prob2_Ritz_Galerk_x_N20');
prob2_y_20 = aysml_read('../dat_dir/prob2_Ritz_Galerk_y_N20');
prob2_v_20 = aysml_read('../dat_dir/prob2_Ritz_Galerk_v_N20');
prob2_w_20 = aysml_read('../dat_dir/prob2_Ritz_Galerk_w_N20');
prob2_u_20 = aysml_read('../dat_dir/prob2_Ritz_Galerk_u_N20');

v_map = @(x, y) x*( 1 - 0.5* (y*y) ).^(0.5) ;
w_map = @(x, y) y*( 1 - 0.5* (x*x) ).^(0.5) ;
u_sol_xy_fun = @(x, y) (1 - v_map(x, y)^2 - w_map(x, y)^2)*exp(- v_map(x, y) );

figure(fig1.Number);
surf(prob2_x_20, prob2_y_20, prob2_u_20)

figure(fig2.Number);
surf(prob2_v_20, prob2_w_20, prob2_u_20)

u_sol = zeros(size(prob2_u_20));

for i=1:size(u_sol, 1)
  for j=1:size(u_sol, 2)
    u_sol(i, j) = u_sol_xy_fun( prob2_x_20(i, j), prob2_y_20(i, j) );
  end
end

figure(fig3.Number)
surf(prob2_x_20, prob2_y_20, u_sol)

figure(fig4.Number)
surf(prob2_v_20, prob2_w_20, u_sol)
