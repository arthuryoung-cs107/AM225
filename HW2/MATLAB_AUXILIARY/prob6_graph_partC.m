fig5 = figure('Name', 'Galaxy Poincare Map, position, problem 6', 'Renderer', 'painters', 'Position', fig_pos(5, :));
xlabel('q1')
ylabel('q3')
hold on

fig6 = figure('Name', 'Galaxy Poincare Map, momentum, problem 6', 'Renderer', 'painters', 'Position', fig_pos(6, :));
xlabel('p1')
ylabel('p2')
zlabel('p3')
view(40, 35)
hold on

prob6_galaxy_sim2 = aysml_read('../dat_dir/prob6_Galaxy_data_sim2');
t_sim2 = prob6_galaxy_sim2(:, 1);
ham_sim2 = prob6_galaxy_sim2(:, 2);
q_sim2 = prob6_galaxy_sim2(:, 3:5);
p_sim2 = prob6_galaxy_sim2(:, 6:8);

q3_positive = q_sim2(:, 3) > 0;
p_sim2_pos = zeros(sum(q3_positive), 3);
p_sim2_neg = zeros(length(q3_positive) - sum(q3_positive), 3);
q_sim2_pos = zeros(sum(q3_positive), 3);
q_sim2_neg = zeros(length(q3_positive) - sum(q3_positive), 3);

pos_count = 1;
neg_count = 1;
for i=1:length(q3_positive)
  if q3_positive(i) == 1
    p_sim2_pos(pos_count, :) = p_sim2(i, :);
    q_sim2_pos(pos_count, :) = q_sim2(i, :);
    pos_count = pos_count + 1;
  else
    p_sim2_neg(neg_count, :) = p_sim2(i, :);
    q_sim2_neg(neg_count, :) = q_sim2(i, :);
    neg_count = neg_count + 1;
  end
end

figure(fig5.Number)
plot(q_sim2_pos(:, 1), q_sim2_pos(:, 3), 'o', 'Color', green4, 'LineWidth', 0.5, 'DisplayName', 'q3 positive')
plot(q_sim2_neg(:, 1), q_sim2_neg(:, 3), 'o', 'Color', red5, 'LineWidth', 0.5, 'DisplayName', 'q3 negative')

figure(fig6.Number)
scatter3( p_sim2_pos(:, 1), p_sim2_pos(:, 2), p_sim2_pos(:, 3), 'Marker', 'o', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', green4, 'LineWidth', 0.5);
scatter3( p_sim2_neg(:, 1), p_sim2_neg(:, 2), p_sim2_neg(:, 3), 'Marker', 'o', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', red5, 'LineWidth', 0.5);
