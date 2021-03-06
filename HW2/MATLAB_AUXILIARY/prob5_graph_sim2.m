fig2 = figure('Name', 'Swarm results, problem 5, sim2', 'Renderer', 'painters', 'Position', fig_pos(1, :));
ylabel('y')
xlabel('x')
xlim([-1, 1])
ylim([-1, 1])
box on
hold on

prob5_part_a_sim2 = aysml_read('../dat_dir/prob5_2DKuramoto_sim2');
% prob5_part_a_sim2 = aysml_read('../dat_dir/prob5_2DKuramoto_sim2_dense_output');


theta_sim2 = zeros(size(prob5_part_a_sim2, 1), dof);
posX_sim2 = zeros(size(prob5_part_a_sim2, 1), dof);
posY_sim2 = zeros(size(prob5_part_a_sim2, 1), dof);
time_sim2 = prob5_part_a_sim2(:, 1);

col_count = 1;
for i=2:3:(size(prob5_part_a_sim2, 2)-2)
  theta_sim2(:, col_count) = prob5_part_a_sim2(:, i);
  posX_sim2(:, col_count) = prob5_part_a_sim2(:, i+1);
  posY_sim2(:, col_count) = prob5_part_a_sim2(:, i+2);
  col_count = col_count + 1;
end

xmax2 = max(max(posX_sim2));
xmin2 = min(min(posX_sim2));
ymax2 = max(max(posY_sim2));
ymin2 = min(min(posY_sim2));

abs2 = max([abs(xmax2), abs(xmin2), abs(ymax2), abs(ymin2)]);

figure(fig2.Number)
for i = 1:length(time_sim2)
  scatter(posX_sim2(i, :)', posY_sim2(i, :)', 30, RGB_theta_gen(theta_sim2(i, :)) );
  ylabel('y')
  xlabel('x')
  xlim([-abs2, abs2])
  ylim([-abs2, abs2])
  box on
  time_sim2(i)
  pause(0.1)
  clf
end
