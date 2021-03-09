fig3 = figure('Name', 'Swarm results, problem 5, sim3', 'Renderer', 'painters', 'Position', fig_pos(1, :));
ylabel('y')
xlabel('x')
xlim([-1, 1])
ylim([-1, 1])
box on
hold on

prob5_part_a_sim3 = aysml_read('../dat_dir/prob5_2DKuramoto_sim3');
% prob5_part_a_sim3 = aysml_read('../dat_dir/prob5_2DKuramoto_sim3_dense_output');

theta_sim3 = zeros(size(prob5_part_a_sim3, 1), dof);
posX_sim3 = zeros(size(prob5_part_a_sim3, 1), dof);
posY_sim3 = zeros(size(prob5_part_a_sim3, 1), dof);
time_sim3 = prob5_part_a_sim3(:, 1);

col_count = 1;
for i=2:3:(size(prob5_part_a_sim3, 2)-2)
  theta_sim3(:, col_count) = prob5_part_a_sim3(:, i);
  posX_sim3(:, col_count) = prob5_part_a_sim3(:, i+1);
  posY_sim3(:, col_count) = prob5_part_a_sim3(:, i+2);
  col_count = col_count + 1;
end

xmax3 = max(max(posX_sim3));
xmin3 = min(min(posX_sim3));
ymax3 = max(max(posY_sim3));
ymin3 = min(min(posY_sim3));

abs3 = max([abs(xmax3), abs(xmin3), abs(ymax3), abs(ymin3)]);

figure(fig3.Number)
for i = 1:length(time_sim3)
  scatter(posX_sim3(i, :)', posY_sim3(i, :)', 30, RGB_theta_gen(theta_sim3(i, :)) );
  ylabel('y')
  xlabel('x')
  xlim([-abs3, abs3])
  ylim([-abs3, abs3])
  box on
  time_sim3(i)
  pause(0.1)
  clf
end
