fig2 = figure('Name', 'Swarm results, problem 5, sim2', 'Renderer', 'painters', 'Position', fig_pos(1, :));
ylabel('y')
xlabel('x')
xlim([-1, 1])
ylim([-1, 1])
box on
hold on

prob5_part_a_sim2 = aysml_read('../dat_dir/prob5_2DKuramoto_sim2');

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

figure(fig2.Number)
for i = 1:length(time_sim2)
  scatter(posX_sim2(i, :)', posY_sim2(i, :)', 30, RGB_theta_gen(theta_sim2(i, :)) );
  ylabel('y')
  xlabel('x')
  xlim([-1.2, 1.2])
  ylim([-1.2, 1.2])
  box on
  pause
  clf
end
