figure1 = figure('Name', 'L 30K', 'Renderer', 'painters', 'Position', fig_pos(1,:));
figure2 = figure('Name', 'L 20K', 'Renderer', 'painters', 'Position', fig_pos(2,:));
figure3 = figure('Name', 'L 10K', 'Renderer', 'painters', 'Position', fig_pos(3,:));
figure5 = figure('Name', 'T 30K', 'Renderer', 'painters', 'Position', fig_pos(5,:));
figure6 = figure('Name', 'T 20K', 'Renderer', 'painters', 'Position', fig_pos(6,:));
figure7 = figure('Name', 'T 10K', 'Renderer', 'painters', 'Position', fig_pos(7,:));

prefix = '../aydat_dir_small_sim';

for i=110:5:110
  for k=1:3
    S_it = (aysml_read([prefix, num2str(2 + k), '/sim_corrupt_', num2str(i), '_S_out']))';
    L_it = (aysml_read([prefix, num2str(2 + k), '/sim_corrupt_', num2str(i), '_L_out']))';
    temp_it = aysml_read([prefix, num2str(2 + k), '/tem', num2str(i)]);
    temp_corrupt = aysml_read([prefix, num2str(2 + k), '/tem_corrupt', num2str(i)]);

    figure(k)
    heatmap(temp_corrupt, 'GridVisible','off');
    Ax1 = gca;
    Ax1.XDisplayLabels = nan(size(Ax1.XDisplayData));
    Ax1.YDisplayLabels = nan(size(Ax1.YDisplayData));

    figure(k+3)
    heatmap(S_it, 'GridVisible','off');
    Ax2 = gca;
    Ax2.XDisplayLabels = nan(size(Ax2.XDisplayData));
    Ax2.YDisplayLabels = nan(size(Ax2.YDisplayData));

    fprintf('sim %d, it = %d \n', k+2, i);
  end
  pause;
end
