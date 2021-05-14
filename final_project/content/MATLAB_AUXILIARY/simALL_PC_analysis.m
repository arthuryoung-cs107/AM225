figure1 = figure('Name', 'mode 1 err', 'Renderer', 'painters', 'Position', fig_pos(1,:));
% set(gca,/ 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('frame')
ylabel('l_2 error')
hold on

figure2 = figure('Name', 'mode 2 err', 'Renderer', 'painters', 'Position', fig_pos(2,:));
% set(gca,/ 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('frame')
ylabel('l_2 error')
hold on

figure3 = figure('Name', 'mode 3 err', 'Renderer', 'painters', 'Position', fig_pos(3,:));
% set(gca,/ 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('frame')
ylabel('l_2 error')
hold on

figure4 = figure('Name', 'mode 4 err', 'Renderer', 'painters', 'Position', fig_pos(4,:));
% set(gca,/ 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('frame')
ylabel('l_2 error')
hold on


template = aysml_read(['../aydat_dir_small_sim3/tem0']);
[U, S, V] = svd(template, 'econ');

U_L_z = zeros(size(U));
S_L_z = zeros(size(S));
V_L_z = zeros(size(V));
U_t_z = zeros(size(U));
S_t_z = zeros(size(S));
V_t_z = zeros(size(V));
U_n_z = zeros(size(U));
S_n_z = zeros(size(S));
V_n_z = zeros(size(V));

err_L = zeros(3, 501, 4);
err_n = zeros(3, 501, 4);

prefix = '../aydat_dir_small_sim';

for k=1:3
  for i=0:1:500
    L_it = (aysml_read([prefix, num2str(2 + k), '/sim_corrupt_', num2str(i), '_L_out']))';
    temp_it = aysml_read([prefix, num2str(2 + k), '/tem', num2str(i)]);
    temp_corrupt = aysml_read([prefix, num2str(2 + k), '/tem_corrupt', num2str(i)]);

    [U_L, S_L, V_L] = svd(L_it, 'econ');
    [U_t, S_t, V_t] = svd(temp_it, 'econ');
    [U_n, S_n, V_n] = svd(temp_corrupt, 'econ');

    for j=1:4
      U_L_z(:, j) = U_L(:, j);
      U_t_z(:, j) = U_t(:, j);
      U_n_z(:, j) = U_n(:, j);
      S_L_z(j, j) = S_L(j, j);
      S_t_z(j, j) = S_t(j, j);
      S_n_z(j, j) = S_n(j, j);
      V_L_z(:, j) = V_L(:, j);
      V_t_z(:, j) = V_t(:, j);
      V_n_z(:, j) = V_n(:, j);

      X_L = U_L_z*S_L_z*(V_L_z');
      X_t = U_t_z*S_t_z*(V_t_z');
      X_n = U_n_z*S_n_z*(V_n_z');

      err_L(k, (i+1), j) = norm(X_L-X_t, 'fro')/(norm(X_t, 'fro'));
      err_n(k, (i+1), j) = norm(X_n-X_t, 'fro')/(norm(X_t, 'fro'));

      U_L_z(:, j) = 0;
      U_t_z(:, j) = 0;
      U_n_z(:, j) = 0;
      S_L_z(j, j) = 0;
      S_t_z(j, j) = 0;
      S_n_z(j, j) = 0;
      V_L_z(:, j) = 0;
      V_t_z(:, j) = 0;
      V_n_z(:, j) = 0;
    end
  end
end

figure(figure1.Number)
plot(0:1:500, err_L(3, :, 1), ' - ', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', '\sigma = 10')
plot(0:1:500, err_L(2, :, 1), ' - ', 'Color', blue5, 'LineWidth', 1.5, 'DisplayName', '\sigma = 20')
plot(0:1:500, err_L(1, :, 1), ' - ', 'Color', green4, 'LineWidth', 1.5, 'DisplayName', '\sigma = 30')
plot(0:1:500, err_n(3, :, 1), ' - ', 'Color', red4, 'LineWidth', 1.5, 'DisplayName', '\sigma = 10 noised')
plot(0:1:500, err_n(2, :, 1), ' - ', 'Color', blue3, 'LineWidth', 1.5, 'DisplayName', '\sigma = 20 noised')
plot(0:1:500, err_n(1, :, 1), ' - ', 'Color', green1, 'LineWidth', 1.5, 'DisplayName', '\sigma = 30 noised')
legend('Show', 'Location', 'SouthWest')

figure(figure2.Number)
plot(0:1:500, err_L(3, :, 2), ' - ', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', '\sigma = 10')
plot(0:1:500, err_L(2, :, 2), ' - ', 'Color', blue5, 'LineWidth', 1.5, 'DisplayName', '\sigma = 20')
plot(0:1:500, err_L(1, :, 2), ' - ', 'Color', green4, 'LineWidth', 1.5, 'DisplayName', '\sigma = 30')
plot(0:1:500, err_n(3, :, 2), ' - ', 'Color', red4, 'LineWidth', 1.5, 'DisplayName', '\sigma = 10 noised')
plot(0:1:500, err_n(2, :, 2), ' - ', 'Color', blue3, 'LineWidth', 1.5, 'DisplayName', '\sigma = 20 noised')
plot(0:1:500, err_n(1, :, 2), ' - ', 'Color', green1, 'LineWidth', 1.5, 'DisplayName', '\sigma = 30 noised')
legend('Show', 'Location', 'SouthWest')

figure(figure3.Number)
plot(0:1:500, err_L(3, :, 3), ' - ', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', '\sigma = 10')
plot(0:1:500, err_L(2, :, 3), ' - ', 'Color', blue5, 'LineWidth', 1.5, 'DisplayName', '\sigma = 20')
plot(0:1:500, err_L(1, :, 3), ' - ', 'Color', green4, 'LineWidth', 1.5, 'DisplayName', '\sigma = 30')
plot(0:1:500, err_n(3, :, 3), ' - ', 'Color', red4, 'LineWidth', 1.5, 'DisplayName', '\sigma = 10 noised')
plot(0:1:500, err_n(2, :, 3), ' - ', 'Color', blue3, 'LineWidth', 1.5, 'DisplayName', '\sigma = 20 noised')
plot(0:1:500, err_n(1, :, 3), ' - ', 'Color', green1, 'LineWidth', 1.5, 'DisplayName', '\sigma = 30 noised')
legend('Show', 'Location', 'SouthWest')

figure(figure4.Number)
plot(0:1:500, err_L(3, :, 4), ' - ', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', '\sigma = 10')
plot(0:1:500, err_L(2, :, 4), ' - ', 'Color', blue5, 'LineWidth', 1.5, 'DisplayName', '\sigma = 20')
plot(0:1:500, err_L(1, :, 4), ' - ', 'Color', green4, 'LineWidth', 1.5, 'DisplayName', '\sigma = 30')
plot(0:1:500, err_n(3, :, 4), ' - ', 'Color', red4, 'LineWidth', 1.5, 'DisplayName', '\sigma = 10 noised')
plot(0:1:500, err_n(2, :, 4), ' - ', 'Color', blue3, 'LineWidth', 1.5, 'DisplayName', '\sigma = 20 noised')
plot(0:1:500, err_n(1, :, 4), ' - ', 'Color', green1, 'LineWidth', 1.5, 'DisplayName', '\sigma = 30 noised')
legend('Show', 'Location', 'SouthWest')
