figure1 = figure('Name', 'S result', 'Renderer', 'painters', 'Position', fig_pos(1,:));
figure2 = figure('Name', 'L result', 'Renderer', 'painters', 'Position', fig_pos(2,:));
figure3 = figure('Name', 'full temp', 'Renderer', 'painters', 'Position', fig_pos(3,:));

prefix = '../aydat_dir_small_sim3/sim1_';

for i=0:1:500

  S_it = (aysml_read([prefix, num2str(i), '_S_out']))';
  L_it = (aysml_read([prefix, num2str(i), '_L_out']))';
  temp_it = aysml_read(['../aydat_dir_small_sim3/tem', num2str(i)]);

  figure(figure1.Number)
  heatmap(S_it)
  figure(figure2.Number)
  heatmap(L_it)
  figure(figure3.Number)
  heatmap(temp_it)

  % fprintf('convergent rank: %d out of %d, error: %f \n', rank(L_out), rank(X_true), norm((L_out-X_true), 'fro')/(norm(X_true, 'fro')) )


  pause(1);
end
