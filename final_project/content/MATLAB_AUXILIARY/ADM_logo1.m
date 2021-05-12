figure1 = figure('Name', 'S result', 'Renderer', 'painters', 'Position', fig_pos(1,:));
figure2 = figure('Name', 'L result', 'Renderer', 'painters', 'Position', fig_pos(2,:));
figure3 = figure('Name', 'pure noise', 'Renderer', 'painters', 'Position', fig_pos(3,:));
figure4 = figure('Name', 'X true', 'Renderer', 'painters', 'Position', fig_pos(4,:));
figure5 = figure('Name', 'X noised', 'Renderer', 'painters', 'Position', fig_pos(7,:));
prefix = '../dat_dir/test2';
% figure6 = figure('Name', 'pure noise, heatmap', 'Renderer', 'painters', 'Position', fig_pos(6,:));


S_out = aysml_read([prefix, '_S_out']);
L_out = aysml_read([prefix, '_L_out']);
pure_noise = aysml_read('../dat_dir/pure_noise');
X_true = aysml_read('../dat_dir/X_true');
X_noised = aysml_read('../dat_dir/X_noised');

S_out = S_out + mean(X_true(:))*ones(size(X_true)); % for improves comparison

S_check = double_to_image2(S_out);
L_check = double_to_image2(L_out);
pure_noise_check = double_to_image2(pure_noise);
X_true_check = double_to_image2(X_true);
X_noised_check = double_to_image2(X_noised);



figure(figure1.Number)
imshow(S_check)
figure(figure2.Number)
imshow(L_check)
figure(figure3.Number)
imshow(pure_noise_check);
figure(figure4.Number)
imshow(X_true_check);
figure(figure5.Number)
imshow(X_noised_check);
% figure(figure6.Number)
% heatmap(pure_noise);

fprintf('convergent rank: %d out of %d, error: %f \n', rank(L_out), rank(X_true), norm((L_out-X_true), 'fro')/(norm(X_true, 'fro')) )
