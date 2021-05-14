figure1 = figure('Name', 'S result', 'Renderer', 'painters', 'Position', fig_pos(1,:));
figure2 = figure('Name', 'L result', 'Renderer', 'painters', 'Position', fig_pos(2,:));
figure3 = figure('Name', 'full temp', 'Renderer', 'painters', 'Position', fig_pos(3,:));
figure4 = figure('Name', 'full temp, corrupted', 'Renderer', 'painters', 'Position', fig_pos(4,:));

prefix = '../aydat_dir_small_sim4/sim_corrupt_';

for i=0:5:500

  S_it = (aysml_read([prefix, num2str(i), '_S_out']))';
  L_it = (aysml_read([prefix, num2str(i), '_L_out']))';
  temp_it = aysml_read(['../aydat_dir_small_sim4/tem', num2str(i)]);
  temp_it_corrupt = aysml_read(['../aydat_dir_small_sim4/tem_corrupt', num2str(i)]);

  figure(figure1.Number)
  heatmap(S_it, 'GridVisible','off');
  Ax1 = gca;
  Ax1.XDisplayLabels = nan(size(Ax1.XDisplayData));
  Ax1.YDisplayLabels = nan(size(Ax1.YDisplayData));

  figure(figure2.Number)
  heatmap(L_it, 'GridVisible','off');
  Ax2 = gca;
  Ax2.XDisplayLabels = nan(size(Ax2.XDisplayData));
  Ax2.YDisplayLabels = nan(size(Ax2.YDisplayData));

  figure(figure3.Number)
  heatmap(temp_it, 'GridVisible','off')
  Ax3 = gca;
  Ax3.XDisplayLabels = nan(size(Ax3.XDisplayData));
  Ax3.YDisplayLabels = nan(size(Ax3.YDisplayData));

  figure(figure4.Number)
  heatmap(temp_it_corrupt, 'GridVisible','off')
  Ax4 = gca;
  Ax4.XDisplayLabels = nan(size(Ax4.XDisplayData));
  Ax4.YDisplayLabels = nan(size(Ax4.YDisplayData));

  % fprintf('convergent rank: %d out of %d, error: %f \n', rank(L_out), rank(X_true), norm((L_out-X_true), 'fro')/(norm(X_true, 'fro')) )


  pause(0.5);
end
