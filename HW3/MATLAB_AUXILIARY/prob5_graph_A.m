fig1 = figure('Name', 'FFT on Poisson, full grid', 'Renderer', 'painters', 'Position', fig_pos(1, :));
xlabel('x')
ylabel('y')
zlabel('v (Schur grid)')
view(40, 35)
hold on

fig2 = figure('Name', 'FFT on Poisson, full grid, source', 'Renderer', 'painters', 'Position', fig_pos(2, :));
xlabel('x')
ylabel('y')
zlabel('f')
view(40, 35)
hold on

fig3 = figure('Name', 'FFT on Poisson, full grid', 'Renderer', 'painters', 'Position', fig_pos(3, :));
xlabel('x')
ylabel('y')
zlabel('v (Schur grid)')
view(40, 35)
hold on

fig4 = figure('Name', 'FFT on Poisson, full grid, source', 'Renderer', 'painters', 'Position', fig_pos(4, :));
xlabel('x')
ylabel('y')
zlabel('f')
view(40, 35)
hold on

fig5 = figure('Name', 'FFT on Poisson, difference', 'Renderer', 'painters', 'Position', fig_pos(5, :));
xlabel('x')
ylabel('y')
zlabel('|GridA - GridB|')
view(40, 35)
hold on

fig6 = figure('Name', 'FFT on Poisson, conjcheck', 'Renderer', 'painters', 'Position', fig_pos(6, :));
xlabel('x')
ylabel('y')
zlabel('v (glue)')
view(40, 35)
hold on


grid_A_conv = aysml_read('../dat_dir/prob5_Poisson_fullgrid');
grid_A_source = aysml_read('../dat_dir/prob5_Poisson_fullgrid_source');

grid_B_conv = aysml_read('../dat_dir/prob5_Poisson_Schur');
grid_B_source = aysml_read('../dat_dir/prob5_Poisson_Schur_source');

grid_B_conjcheck = aysml_read('../dat_dir/prob5_Poisson_Schur_conjcheck');

figure(fig1.Number);
surf(grid_A_conv);

figure(fig2.Number);
surf(grid_A_source);

figure(fig3.Number);
surf(grid_B_conv);

figure(fig4.Number);
surf(grid_B_source);

err_mat = grid_A_conv-grid_B_conv;
figure(fig5.Number);
surf(abs(err_mat));

figure(fig6.Number);
surf(grid_B_conjcheck);
