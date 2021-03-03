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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig1 = figure('Name', 'Brusselator results, problem 2', 'Renderer', 'painters', 'Position', fig_pos(1, :));
ylabel('efficiency')
% xlabel('n')
box on
hold on

prob2_part_a_3 = aysml_read('../dat_dir/prob2_Bruss_results_lambda-3');
prob2_part_a_4 = aysml_read('../dat_dir/prob2_Bruss_results_lambda-4');
prob2_part_a_5 = aysml_read('../dat_dir/prob2_Bruss_results_lambda-5');
prob2_part_a_6 = aysml_read('../dat_dir/prob2_Bruss_results_lambda-6');
prob2_part_a_7 = aysml_read('../dat_dir/prob2_Bruss_results_lambda-7');
prob2_part_a_8 = aysml_read('../dat_dir/prob2_Bruss_results_lambda-8');
prob2_part_a_9 = aysml_read('../dat_dir/prob2_Bruss_results_lambda-9');
prob2_part_a_10 = aysml_read('../dat_dir/prob2_Bruss_results_lambda-10');
prob2_part_a_11 = aysml_read('../dat_dir/prob2_Bruss_results_lambda-11');
prob2_part_a_12 = aysml_read('../dat_dir/prob2_Bruss_results_lambda-12');
prob2_part_a_13 = aysml_read('../dat_dir/prob2_Bruss_results_lambda-13');

figure(fig1.Number)
plot(prob2_part_a_3(:, 1), prob2_part_a_3(:, 2), '- o', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_3(:, 1), prob2_part_a_3(:, 3), '- s', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_4(:, 1), prob2_part_a_4(:, 2), '- o', 'Color', pink(1, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_4(:, 1), prob2_part_a_4(:, 3), '- s', 'Color', pink(1, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_5(:, 1), prob2_part_a_5(:, 2), '- o', 'Color', pink(2, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_5(:, 1), prob2_part_a_5(:, 3), '- s', 'Color', pink(2, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_6(:, 1), prob2_part_a_6(:, 2), '- o', 'Color', pink(3, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_6(:, 1), prob2_part_a_6(:, 3), '- s', 'Color', pink(3, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_7(:, 1), prob2_part_a_7(:, 2), '- o', 'Color', pink(4, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_7(:, 1), prob2_part_a_7(:, 3), '- s', 'Color', pink(4, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_8(:, 1), prob2_part_a_8(:, 2), '- o', 'Color', pink(5, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_8(:, 1), prob2_part_a_8(:, 3), '- s', 'Color', pink(5, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_9(:, 1), prob2_part_a_9(:, 2), '- o', 'Color', pink(6, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_9(:, 1), prob2_part_a_9(:, 3), '- s', 'Color', pink(6, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_10(:, 1), prob2_part_a_10(:, 2), '- o', 'Color', pink(7, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_10(:, 1), prob2_part_a_10(:, 3), '- s', 'Color', pink(7, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_11(:, 1), prob2_part_a_11(:, 2), '- o', 'Color', pink(8, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_11(:, 1), prob2_part_a_11(:, 3), '- s', 'Color', pink(8, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_12(:, 1), prob2_part_a_12(:, 2), '- o', 'Color', pink(9, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_12(:, 1), prob2_part_a_12(:, 3), '- s', 'Color', pink(9, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_13(:, 1), prob2_part_a_13(:, 2), '- o', 'Color', pink(10, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(prob2_part_a_13(:, 1), prob2_part_a_13(:, 3), '- s', 'Color', pink(10, :), 'LineWidth', 1.5, 'DisplayName', '2 threads')
