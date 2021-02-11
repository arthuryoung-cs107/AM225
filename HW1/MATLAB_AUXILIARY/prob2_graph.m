clear
close all

fig_pos = [10, 455, 350, 350;
    360, 455, 350, 350; 710, 455, 350, 350; 1060, 455, 350, 350; 10, 30, 350, 350; 360, 30, 350, 350; 710, 30, 350, 350; 1060, 30, 350, 350];

green1 = [88/250, 214/250, 141/250];
green2 = [40/250, 180/250, 99/250];
green3 = [34/250, 153/250, 84/250];
green4 = [25/250, 111/250, 61/250];
green5 = [0, 1, 0];
green = [green1; green2; green3; green4; green5];

blue1 = [120/250, 150/250, 250/250];
blue2 = [52/250, 152/250, 219/250];
blue3 = [39/250, 97/250, 141/250];
blue4 = [10/250, 50/250, 150/250];
blue5 = [0, 0, 1];
blue = [blue1; blue2; blue3; blue4; blue5];

red1 = [236/250, 112/250, 99/250];
red2 = [192/250, 57/250, 43/250];
red3 = [146/250, 43/250, 33/250];
red4 = [100/250, 30/250 , 22/250];
red5 = [1, 0 , 0];
red = [red1; red2; red3; red4; red5];

pink = [255 0 104 ; 243 0 112; 230 0 119 ; 216  0 125; 200 0 131; 183 0 136; 165 0 140; 146 0 143; 125 7 145; 102 20 146]*(255)^(-1);

view_mat = [45, 45; 1, 0; 0, 90; 90, 0 ; 45, 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dims = dlmread('m_n.dat');
m = dims(1);
n = dims(2);

gradients_analytical = grad_eval_analytical(points);

fig1 = figure('Name', 'Jet space, view 1', 'Renderer', 'painters', 'Position', fig_pos(1, :));
view(view_mat(1, :));
ylabel('x')
xlabel('t')
xlim([0 6])
ylim([0 9])
hold on

sol_it = 1;

figure(fig1.Number)
scatter3(points(:, 1), points(:, 2), points(:, 3), ' o','LineWidth', 1, 'CData', points(:, 3))

figure(fig2.Number)
scatter3(points_noise(:, 1), points_noise(:, 2), points_noise(:, 3), ' o','LineWidth', 1, 'CData', points(:, 3))

% un-noised G matrix
G_true_visible = dlmread('G_true_visible.dat');
G_true = (fread(fopen('./svdTRUE_dir/G_true.dat'), [n, m], 'float64=>float64'))';
