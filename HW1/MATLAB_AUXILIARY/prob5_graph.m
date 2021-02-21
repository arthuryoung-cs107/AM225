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

% figure
% mesh(aysml_read('../interim_matrix_check'));

fig1 = figure('Name', 'Jet space, view 1', 'Renderer', 'painters', 'Position', fig_pos(1, :));
ylabel('Wall time (seconds)')
xlabel('N')
box on
hold on

fig2 = figure('Name', 'Jet space, view 1', 'Renderer', 'painters', 'Position', fig_pos(2, :));
ylabel('S')
xlabel('time')
box on
hold on


thread1_results = aysml_read('../dat_dir/prob5_threads1_runtime');
thread2_results = aysml_read('../dat_dir/prob5_threads2_runtime');
thread3_results = aysml_read('../dat_dir/prob5_threads3_runtime');
thread4_results = aysml_read('../dat_dir/prob5_threads4_runtime');
thread8_results = aysml_read('../dat_dir/prob5_threads8_runtime');

s_vs_t_T1 = aysml_read('../dat_dir/prob5_threads8_s_vs_t_T1');
s_vs_t_Tid = aysml_read('../dat_dir/prob5_threads8_s_vs_t_Tid');


figure(fig1.Number)
plot(thread1_results(:, 2), thread1_results(:, 3), '- o', 'Color', [0 0 0], 'LineWidth', 1.5, 'DisplayName', '1 thread')
plot(thread2_results(:, 2), thread2_results(:, 3), '- o', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(thread3_results(:, 2), thread3_results(:, 3), '- o', 'Color', green4, 'LineWidth', 1.5, 'DisplayName', '3 threads')
plot(thread4_results(:, 2), thread4_results(:, 3), '- o', 'Color', blue5, 'LineWidth', 1.5, 'DisplayName', '4 threads')
% plot(thread8_results(:, 2), thread8_results(:, 3), '- o', 'Color', pink(1, :), 'LineWidth', 1.5, 'DisplayName', '8 threads')

figure(fig2.Number)
plot(s_vs_t_T1(:, 1), s_vs_t_T1(:, 2), '- .', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', 'T original')
plot(s_vs_t_Tid(:, 1), s_vs_t_Tid(:, 2), '- .', 'Color', blue5, 'LineWidth', 1.5, 'DisplayName', 'T identity')

legend('Show', 'Location', 'NorthEast')
