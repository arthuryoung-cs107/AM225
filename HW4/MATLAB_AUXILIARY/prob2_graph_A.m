fig1 = figure('Name', 'FFT on Poisson, full grid', 'Renderer', 'painters', 'Position', fig_pos(1, :));
xlabel('x')
ylabel('y')
zlabel('v')
view(40, 35)
hold on

fig2 = figure('Name', 'FFT on Poisson, full grid, source', 'Renderer', 'painters', 'Position', fig_pos(2, :));
xlabel('x')
ylabel('y')
zlabel('f')
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

grid_A_conv = aysml_read('../dat_dir/prob2_square_fft');
grid_A_source = aysml_read('../dat_dir/prob2_square_fft_source');

figure(fig1.Number);
surf(grid_A_conv);

figure(fig2.Number);
surf(grid_A_source);

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
