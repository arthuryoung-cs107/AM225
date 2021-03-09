fig3 = figure('Name', 'Galaxy sim1, problem 6', 'Renderer', 'painters', 'Position', fig_pos(3, :));
xlabel('q1')
ylabel('q2')
zlabel('q3')
view(40, 35)
hold on

fig4 = figure('Name', 'Galaxy Hamiltonian, problem 6', 'Renderer', 'painters', 'Position', fig_pos(4, :));
ylabel('H')
xlabel('time (dimensionless)')
box on
hold on


prob6_galaxy_sim1 = aysml_read('../dat_dir/prob6_Galaxy_data_sim1');
t_sim1 = prob6_galaxy_sim1(:, 1);
ham_sim1 = prob6_galaxy_sim1(:, 2);
q_sim1 = prob6_galaxy_sim1(:, 3:5);
p_sim1 = prob6_galaxy_sim1(:, 6:8);
deltaT = length(t_sim1)/100;

figure(fig4.Number)
plot(t_sim1, ham_sim1, '- ', 'Color', blue5, 'LineWidth', 1.5, 'DisplayName', 'Hamiltonian value')

pause()


xmax1 = max(max(q_sim1));
xmin1 = min(min(q_sim1));

abs1 = max([abs(xmax1), abs(xmax1)]);

figure(fig3.Number)
for i=1:2:length(t_sim1)
  % scatter3( q_sim1( 1, 1), q_sim1(1, 2), q_sim1(1, 3), 500,'Marker', 'x', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 2);
  scatter3( q_sim1( i, 1), q_sim1(i, 2), q_sim1(i, 3), 'Marker', 'o', 'MarkerEdgeColor', green4, 'MarkerFaceColor', blue2, 'LineWidth', 0.5);
  % scatter3( q_sim1( i + deltaT-1, 1), q_sim1(i + deltaT-1 , 2), q_sim1(i + deltaT-1 , 3), 500, 'filled', 'Marker', '*', 'MarkerEdgeColor', red5, 'LineWidth', 2);
  xlim([-abs1, abs1])
  ylim([-abs1, abs1])
  zlim([-abs1, abs1])
  xlabel('q1')
  ylabel('q2')
  zlabel('q3')
  view(40, 35)
  t_sim1(i)
  pause(0.001)
  clf
end
