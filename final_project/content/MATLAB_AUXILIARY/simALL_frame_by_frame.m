figure1 = figure('Name', 'S result', 'Renderer', 'painters', 'Position', fig_pos(1,:));
figure2 = figure('Name', 'L result', 'Renderer', 'painters', 'Position', fig_pos(2,:));
figure3 = figure('Name', 'full temp', 'Renderer', 'painters', 'Position', fig_pos(3,:));
figure4 = figure('Name', 'full temp, corrupted', 'Renderer', 'painters', 'Position', fig_pos(4,:));

prefix = '../aydat_dir_small_sim';

for k=1:3
  for i=50:25:300
    S_it = (aysml_read([prefix, num2str(2 + k), '/sim_corrupt_', num2str(i), '_S_out']))';
    L_it = (aysml_read([prefix, num2str(2 + k), '/sim_corrupt_', num2str(i), '_L_out']))';
    temp_it = aysml_read([prefix, num2str(2 + k), '/tem', num2str(i)]);
    temp_corrupt = aysml_read([prefix, num2str(2 + k), '/tem_corrupt', num2str(i)]);

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
    heatmap(temp_corrupt, 'GridVisible','off')
    Ax4 = gca;
    Ax4.XDisplayLabels = nan(size(Ax4.XDisplayData));
    Ax4.YDisplayLabels = nan(size(Ax4.YDisplayData));

    fprintf('sim %d, it = %d \n', k+2, i);

    % pause;
    pause(1);
  end
end
