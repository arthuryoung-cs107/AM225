fig2 = figure('Name', 'Swarm results, problem 5, sim56', 'Renderer', 'painters', 'Position', fig_pos(1, :));
ylabel('y')
xlabel('x')
xlim([-1, 1])
ylim([-1, 1])
box on
hold on

prob5_part_a_sim56 = aysml_read('../dat_dir/prob5_2DKuramotoSuper_sim5-6');


theta_sim56 = zeros(size(prob5_part_a_sim56, 1), dof);
posX_sim56 = zeros(size(prob5_part_a_sim56, 1), dof);
posY_sim56 = zeros(size(prob5_part_a_sim56, 1), dof);
time_sim56 = prob5_part_a_sim56(:, 1);

col_count = 1;
for i=2:3:(size(prob5_part_a_sim56, 2)-2)
  theta_sim56(:, col_count) = prob5_part_a_sim56(:, i);
  posX_sim56(:, col_count) = prob5_part_a_sim56(:, i+1);
  posY_sim56(:, col_count) = prob5_part_a_sim56(:, i+2);
  col_count = col_count + 1;
end

xmax56 = max(max(posX_sim56));
xmin56 = min(min(posX_sim56));
ymax56 = max(max(posY_sim56));
ymin56 = min(min(posY_sim56));

abs56 = max([abs(xmax56), abs(xmin56), abs(ymax56), abs(ymin56)]);

figure(fig2.Number)
for i = 1:length(time_sim56)
  scatter(posX_sim56(i, :)', posY_sim56(i, :)', 30, RGB_theta_gen(theta_sim56(i, :)) );
  ylabel('y')
  xlabel('x')
  xlim([-abs56, abs56])
  ylim([-abs56, abs56])
  box on
  time_sim56(i)
  pause(0.01)
  clf
end
