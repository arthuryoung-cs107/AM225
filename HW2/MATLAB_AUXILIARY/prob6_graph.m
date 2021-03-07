clear
close all

fig_pos = [10, 455, 350, 350;
    360, 455, 350, 350; 710, 455, 350, 350; 1060, 455, 350, 350; 10, 30, 350, 350; 360, 30, 350, 350; 710, 30, 350, 350; 1060, 30, 350, 350];

grey1 = [178/250, 186/250, 187/250];
grey2 = [131/250, 145/250, 146/250];
grey3 = [97/250, 106/250, 107/250];
grey4 = [66/250, 73/250, 73/250];
grey5 = [20/100 20/100 20/100];
grey = [grey1; grey2; grey3; grey4; grey5];

purple1 = [102/250, 0/250, 102/250];
purple2 = [153/250, 0/250, 153/250];
purple3 = [204/250, 0/250, 204/250];
purple4 = [250/250, 0/250, 250/250];
purple5 = [250/250, 50/250, 250/250];
purple = [purple1; purple2; purple3; purple4; purple5];

orange1 = [255/255 90/255 0];
orange2 = [255/255 123/255 0];
orange3 = [255/255 165/255 0];
orange4 = [255/255 208/255 0];
orange5 = [255/255 229/255 0];
orange = [orange1; orange2; orange3; orange4; orange5];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig1 = figure('Name', 'Brusselator results, problem 6', 'Renderer', 'painters', 'Position', fig_pos(1, :));
ylabel('y')
xlabel('t')
box on
hold on

fig2 = figure('Name', 'Brusselator Convergence, problem 6', 'Renderer', 'painters', 'Position', fig_pos(2, :));
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold on
ylabel('l2 error')
xlabel('delT')
box on


prob6_parta_Brus_1e2 = aysml_read('../dat_dir/prob6_BrussGeng_results_steps_1e2');
prob6_parta_Brus_5e2 = aysml_read('../dat_dir/prob6_BrussGeng_results_steps_5e2');
prob6_parta_Brus_1e3 = aysml_read('../dat_dir/prob6_BrussGeng_results_steps_1e3');
prob6_parta_Brus_5e3 = aysml_read('../dat_dir/prob6_BrussGeng_results_steps_5e3');
prob6_parta_Brus_1e4 = aysml_read('../dat_dir/prob6_BrussGeng_results_steps_1e4');
prob6_parta_Brus_5e4 = aysml_read('../dat_dir/prob6_BrussGeng_results_steps_5e4');
prob6_parta_Brus_1e5 = aysml_read('../dat_dir/prob6_BrussGeng_results_steps_1e5');
prob6_parta_Brus_5e5 = aysml_read('../dat_dir/prob6_BrussGeng_results_steps_5e5');
prob6_parta_Brus_1e6 = aysml_read('../dat_dir/prob6_BrussGeng_results_steps_1e6');
prob6_parta_Brus_5e6 = aysml_read('../dat_dir/prob6_BrussGeng_results_steps_5e6');
prob6_parta_Brus_1e7 = aysml_read('../dat_dir/prob6_BrussGeng_results_steps_1e7');
% prob6_parta_Brus_5e7 = aysml_read('../dat_dir/prob6_BrussGeng_results_steps_5e7');
% prob6_parta_Brus_1e8 = aysml_read('../dat_dir/prob6_BrussGeng_results_steps_1e8');

prob2_part_a_15 = aysml_read('../dat_dir/prob2_Bruss_results_lambda-15');

% precision_vec = zeros(13, 1);
precision_vec = zeros(11, 1);

precision_vec(1) = precision_func(prob6_parta_Brus_1e2, prob2_part_a_15);
precision_vec(2) = precision_func(prob6_parta_Brus_5e2, prob2_part_a_15);
precision_vec(3) = precision_func(prob6_parta_Brus_1e3, prob2_part_a_15);
precision_vec(4) = precision_func(prob6_parta_Brus_5e3, prob2_part_a_15);
precision_vec(5) = precision_func(prob6_parta_Brus_1e4, prob2_part_a_15);
precision_vec(6) = precision_func(prob6_parta_Brus_5e4, prob2_part_a_15);
precision_vec(7) = precision_func(prob6_parta_Brus_1e5, prob2_part_a_15);
precision_vec(8) = precision_func(prob6_parta_Brus_5e5, prob2_part_a_15);
precision_vec(9) = precision_func(prob6_parta_Brus_1e6, prob2_part_a_15);
precision_vec(10) = precision_func(prob6_parta_Brus_5e6, prob2_part_a_15);
precision_vec(11) = precision_func(prob6_parta_Brus_1e7, prob2_part_a_15);
% precision_vec(12) = precision_func(prob6_parta_Brus_5e7, prob2_part_a_15);
% precision_vec(13) = precision_func(prob6_parta_Brus_1e8, prob2_part_a_15);

delT_vec = zeros(size(precision_vec));

delT_vec(1)= 1e2;
delT_vec(2)= 5e2;
delT_vec(3)= 1e3;*.dat
delT_vec(4)= 5e3;
delT_vec(5)= 1e4;
delT_vec(6)= 5e4;
delT_vec(7)= 1e5;
delT_vec(8) = 5e5;
delT_vec(9) = 1e6;
delT_vec(10) = 5e6;
delT_vec(11) = 1e7;
% delT_vec(12) = 5e7;
% delT_vec(13) = 1e8;

delT_vec = (20.0)./( delT_vec );

figure(fig1.Number)
plot(prob6_parta_Brus_1e7(:, 1), prob6_parta_Brus_1e7(:, 2), '- ', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', 'y1: delT = 1e-1')
plot(prob6_parta_Brus_1e7(:, 1), prob6_parta_Brus_1e7(:, 3), '- ', 'Color', blue5, 'LineWidth', 1.5, 'DisplayName', 'y2: lamba = 1e-1')

figure(fig2.Number)
plot(delT_vec, precision_vec, '- ', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', 'convergence')
