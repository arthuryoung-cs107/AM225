close all
clear

fig_pos_old = [10, 455, 350, 350; 360, 455, 350, 350; 710, 455, 350, 350; 1060, 455, 350, 350; 10, 30, 350, 350; 360, 30, 350, 350; 710, 30, 350, 350; 1060, 30, 350, 350];
fig_pos = [1 498 173 130; 175 498 173 130; 349 498 173 130; 523 498 173 130; 697 498 173 130; 871 498 173 130; 1045 498 173 130];

figure1 = figure('Name', 'Imported Image', 'Renderer', 'painters', 'Position', fig_pos(1,:));
figure2 = figure('Name', 'Grey scale check', 'Renderer', 'painters', 'Position', fig_pos(2,:));
figure3 = figure('Name', 'Truncated Rank', 'Renderer', 'painters', 'Position', fig_pos(3,:));
figure4 = figure('Name', 'Noised Truncated Rank', 'Renderer', 'painters', 'Position', fig_pos(4,:));
figure5 = figure('Name', 'convergence check', 'Renderer', 'painters', 'Position', fig_pos(5,:));
figure6 = figure('Name', 'convergence check2', 'Renderer', 'painters', 'Position', fig_pos(6,:));


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
S_0 = S_00; %% truncate the rank
S_0(r+1:n, r+1:n) = 0;
X0 = U_00* S_0* V_00';
X0_check = double_to_image2(X0);
figure(figure3.Number)
imshow(X0_check)

X0_vec = X0(:);
rand_mat = randi([0, 1], [m, n]);
Omega = find(rand_mat);
M_raw = zeros(length(Omega), 1);
X = zeros(m*n, 1);
count = 1;
for i = Omega'
    M_raw(count) = X0_vec(i);
    X(i) = X0_vec(i);
    count = count + 1;
end
X = reshape(X, [m, n]);

fprintf_matrix(X, '../dat_dir/X00')

[U,S,V] = svd(X, 'econ');
X_check1 = double_to_image2(X);
figure(figure4.Number)
imshow(X_check1)
X = X/S(1, 1);
M = M_raw/S(1, 1);

fprintf_matrix(X, '../dat_dir/X0')

mat_test = [1 2 3 4 ; 5 6 7 8 ; 9 10 11 12];
fprintf_matrix(mat_test, '../dat_dir/test')
