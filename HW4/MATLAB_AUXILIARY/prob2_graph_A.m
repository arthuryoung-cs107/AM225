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

Ritz_galerk_xyvw = aysml_read('../dat_dir/prob2_Ritz_Galerk_xyvw');
Ritz_galerk_u = aysml_read('../dat_dir/prob2_Ritz_Galerk_u');

figure(fig1.Number);
surf(Ritz_galerk_xyvw(1:(size(Ritz_galerk_u, 1)), 2), Ritz_galerk_xyvw(1:(size(Ritz_galerk_u, 1)), 2), Ritz_galerk_u)
%
% figure(fig2.Number);
% surf(xyvw_mat(:, 3), xyvw_mat(:,4), phi_mat)

x_vec = -1:0.1:1;
y_vec = -1:0.1:1;

v_map = @(x, y) x*( 1 - 0.5* (y*y) ).^(0.5) ;
w_map = @(x, y) y*( 1 - 0.5* (x*x) ).^(0.5) ;

u_sol_xy_fun = @(x, y) (1 - v_map(x, y)^2 - w_map(x, y)^2)*exp(- v_map(x, y) );
u_sol = zeros(length(y_vec), length(x_vec) );
v_mat = zeros(length(y_vec), length(x_vec) );
w_mat = zeros(length(y_vec), length(x_vec) );


for i=1:length(y_vec)
  for j=1:length(x_vec)
    u_sol(i, j) = u_sol_xy_fun( x_vec(j), y_vec(i) );
    v_mat(i, j) = v_map( x_vec(j), y_vec(i) );
    w_mat(i, j) = w_map( x_vec(j), y_vec(i) );
  end
end

figure(fig3.Number)
surf(x_vec, y_vec, u_sol);

figure(fig4.Number)
surf(w_mat, v_mat, u_sol);
