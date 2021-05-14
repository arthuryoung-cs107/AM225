figure1 = figure('Name', 'Imported Image', 'Renderer', 'painters', 'Position', fig_pos(1,:));
figure2 = figure('Name', 'Grey scale check', 'Renderer', 'painters', 'Position', fig_pos(2,:));
figure3 = figure('Name', 'Truncated Rank', 'Renderer', 'painters', 'Position', fig_pos(3,:));
figure4 = figure('Name', 'Noised Truncated Rank', 'Renderer', 'painters', 'Position', fig_pos(4,:));
figure5 = figure('Name', 'noised matrix', 'Renderer', 'painters', 'Position', fig_pos(5,:));
figure6 = figure('Name', 'pure noise', 'Renderer', 'painters', 'Position', fig_pos(6,:));


X_rgb = imread('harvard_logo.png');
X_grey_uint8 = rgb2gray(X_rgb);
figure(figure1.Number)
imshow(X_rgb)
X_grey = double(X_grey_uint8);
s_image = svd(X_grey);

X_grey_check = double_to_image2(X_grey);
figure(figure2.Number)
imshow(X_grey_check)

X00 = X_grey;
[U_00,S_00,V_00]=svd(X_grey, 'econ');
[m, n] = size(X00);
r00 = rank(X00);
r = ceil(r00/6);
% r = ceil(r00/2);
% r = r00;
S_0 = S_00; %% truncate the rank
S_0(r+1:n, r+1:n) = 0;
X0 = U_00* S_0* V_00';
X0_check = double_to_image2(X0);
figure(figure3.Number)
imshow(X0_check)

fprintf_matrix(X0, '../logo_dat_dir/X_true')
fprintf_matrix(X0/S_0(1), '../logo_dat_dir/Xnorm_true')

X0_vec = X0(:);

mu = mean(X0_vec);
sigma = std(X0_vec);

rand_mat = (unifrnd(0, 1, [m, n])) > 0.1; % percent sparsity
Omega = find(rand_mat); % find un-noised values
Omega_not = find(rand_mat == 0); % find noised values
bin_rand = 2*(randi([0, 1], size(Omega_not))) - 1;
uni_rand = 255*unifrnd(-0.5, 0.5, size(Omega_not));
M_raw = zeros(length(Omega), 1);
X = zeros(m*n, 1);
pure_noise = mu*ones(m*n, 1);
Xn = X0_vec;
count = 1;
for i = Omega'
    M_raw(count) = X0_vec(i);
    X(i) = X0_vec(i);
    count = count + 1;
end
count = 1;
for i = Omega_not'
  rand_val = 10*sigma*bin_rand(count);
  % rand_val = uni_rand(count);
  Xn(i) = Xn(i) + rand_val;
  % Xn(i) = mu + rand_val;
  pure_noise(i) = rand_val;
  count = count + 1;
end
X = reshape(X, [m, n]);
Xn = reshape(Xn, [m, n]);
pure_noise = reshape(pure_noise, [m, n]);


fprintf_matrix(X, '../logo_dat_dir/X_completion')
fprintf_matrix(Xn, '../logo_dat_dir/X_noised')
fprintf_matrix(pure_noise, '../logo_dat_dir/pure_noise')

[U,S,V] = svd(X, 'econ');
[Un,Sn,Vn] = svd(Xn, 'econ');
X_check1 = double_to_image2(X);
Xn_check1 = double_to_image2(Xn);
pure_noise_check1 = double_to_image2(pure_noise);
figure(figure4.Number)
imshow(X_check1)
figure(figure5.Number)
imshow(Xn_check1)
figure(figure6.Number)
imshow(pure_noise_check1)


X = X/S(1, 1);
Xn = Xn/Sn(1, 1);
M = M_raw/S(1, 1);

fprintf_matrix(X, '../logo_dat_dir/Xnorm_completion')
fprintf_matrix(Xn, '../logo_dat_dir/Xnorm_noised')
