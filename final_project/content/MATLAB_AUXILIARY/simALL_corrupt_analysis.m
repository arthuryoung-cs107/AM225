
figure1 = figure('Name', 'err results', 'Renderer', 'painters', 'Position', fig_pos(1,:));
% set(gca,/ 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('frame')
ylabel('l_2 error')
hold on

figure2 = figure('Name', 'sigma results', 'Renderer', 'painters', 'Position', fig_pos(2,:));
% set(gca,/ 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('frame')
ylabel('\sigma(\chi)')
hold on


prefix3 = '../aydat_dir_small_sim3/sim_corrupt_';
prefix4 = '../aydat_dir_small_sim4/sim_corrupt_';
prefix5 = '../aydat_dir_small_sim5/sim_corrupt_';

err_3 = zeros(1, 501);
err_4 = zeros(1, 501);
err_5 = zeros(1, 501);

err_n_3 = zeros(1, 501);
err_n_4 = zeros(1, 501);
err_n_5 = zeros(1, 501);

std_3 = zeros(1, 501);
std_4 = zeros(1, 501);
std_5 = zeros(1, 501);

for i=0:1:500

  L_it_3 = (aysml_read([prefix3, num2str(i), '_L_out']))';
  temp_it_3 = aysml_read(['../aydat_dir_small_sim3/tem', num2str(i)]);
  temp_it_corrupt_3 = aysml_read(['../aydat_dir_small_sim3/tem_corrupt', num2str(i)]);

  L_it_4 = (aysml_read([prefix4, num2str(i), '_L_out']))';
  temp_it_4 = aysml_read(['../aydat_dir_small_sim4/tem', num2str(i)]);
  temp_it_corrupt_4 = aysml_read(['../aydat_dir_small_sim4/tem_corrupt', num2str(i)]);

  L_it_5 = (aysml_read([prefix5, num2str(i), '_L_out']))';
  temp_it_5 = aysml_read(['../aydat_dir_small_sim5/tem', num2str(i)]);
  temp_it_corrupt_5 = aysml_read(['../aydat_dir_small_sim5/tem_corrupt', num2str(i)]);

  err_3(i + 1) = norm(L_it_3-temp_it_3, 'fro')/norm(temp_it_3, 'fro');
  err_4(i + 1) = norm(L_it_4-temp_it_4, 'fro')/norm(temp_it_4, 'fro');
  err_5(i + 1) = norm(L_it_5-temp_it_5, 'fro')/norm(temp_it_5, 'fro');

  err_n_3(i + 1) = norm(temp_it_corrupt_3-temp_it_3, 'fro')/norm(temp_it_3, 'fro');
  err_n_4(i + 1) = norm(temp_it_corrupt_4-temp_it_4, 'fro')/norm(temp_it_4, 'fro');
  err_n_5(i + 1) = norm(temp_it_corrupt_5-temp_it_5, 'fro')/norm(temp_it_5, 'fro');

  std_3(i + 1) = std(temp_it_3(:));
  std_4(i + 1) = std(temp_it_4(:));
  std_5(i + 1) = std(temp_it_5(:));

end

figure(figure1.Number)
plot(0:1:500, err_5, ' - ', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', '\sigma = 10 ')
plot(0:1:500, err_4, ' - ', 'Color', blue5, 'LineWidth', 1.5, 'DisplayName', '\sigma = 20 ')
plot(0:1:500, err_3, ' - ', 'Color', green4, 'LineWidth', 1.5, 'DisplayName', '\sigma = 30 ')
plot(0:1:500, err_n_5, ' - ', 'Color', red4, 'LineWidth', 1.5, 'DisplayName', '\sigma = 10 noised ')
plot(0:1:500, err_n_4, ' - ', 'Color', blue3, 'LineWidth', 1.5, 'DisplayName', '\sigma = 20 noised ')
plot(0:1:500, err_n_3, ' - ', 'Color', green1, 'LineWidth', 1.5, 'DisplayName', '\sigma = 30 noised ')
legend('Show', 'Location', 'SouthWest')

figure(figure2.Number)
plot(0:1:500, std_5, ' - ', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', '\sigma = 10 ')
plot(0:1:500, std_4, ' - ', 'Color', blue5, 'LineWidth', 1.5, 'DisplayName', '\sigma = 20 ')
plot(0:1:500, std_3, ' - ', 'Color', green4, 'LineWidth', 1.5, 'DisplayName', '\sigma = 30 ')
legend('Show', 'Location', 'SouthWest')
