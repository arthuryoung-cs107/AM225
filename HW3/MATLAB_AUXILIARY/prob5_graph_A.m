fig1 = figure('Name', 'FFT on Poisson, full grid', 'Renderer', 'painters', 'Position', fig_pos(1, :));
xlabel('x')
ylabel('y')
zlabel('v')
view(40, 35)
hold on

fig2 = figure('Name', 'FFT on Poisson, full grid, source', 'Renderer', 'painters', 'Position', fig_pos(2, :));
xlabel('x')
ylabel('y')
zlabel('v')
view(40, 35)
hold on


grid_A_conv = aysml_read('../dat_dir/prob5_Poisson_fullgrid');
grid_A_source = aysml_read('../dat_dir/prob5_Poisson_fullgrid_source');

figure(fig1.Number);
surf(grid_A_conv);

figure(fig2.Number);
surf(grid_A_source);
