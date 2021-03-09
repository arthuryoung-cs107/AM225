fig1 = figure('Name', 'Swarm results, problem 5, sim1', 'Renderer', 'painters', 'Position', fig_pos(1, :));
ylabel('y')
xlabel('x')
xlim([-1, 1])
ylim([-1, 1])
box on
hold on

prob5_part_a_sim1 = aysml_read('../dat_dir/prob5_2DKuramoto_sim1');
% prob5_part_a_sim1 = aysml_read('../dat_dir/prob5_2DKuramoto_sim1_dense_output');

theta_sim1 = zeros(size(prob5_part_a_sim1, 1), dof);
posX_sim1 = zeros(size(prob5_part_a_sim1, 1), dof);
posY_sim1 = zeros(size(prob5_part_a_sim1, 1), dof);
time_sim1 = prob5_part_a_sim1(:, 1);

col_count = 1;
for i=2:3:(size(prob5_part_a_sim1, 2)-2)
  theta_sim1(:, col_count) = prob5_part_a_sim1(:, i);
  posX_sim1(:, col_count) = prob5_part_a_sim1(:, i+1);
  posY_sim1(:, col_count) = prob5_part_a_sim1(:, i+2);
  col_count = col_count + 1;
end

xmax1 = max(max(posX_sim1));
xmin1 = min(min(posX_sim1));
ymax1 = max(max(posY_sim1));
ymin1 = min(min(posY_sim1));

abs1 = max([abs(xmax1), abs(xmin1), abs(ymax1), abs(ymin1)]);

figure(fig1.Number)
for i = 1:length(time_sim1)
  scatter(posX_sim1(i, :)', posY_sim1(i, :)', 30, RGB_theta_gen(theta_sim1(i, :)) );
  ylabel('y')
  xlabel('x')
  xlim([-abs1, abs1])
  ylim([-abs1, abs1])
  box on
  time_sim1(i)
  pause(0.1)
  clf
end
