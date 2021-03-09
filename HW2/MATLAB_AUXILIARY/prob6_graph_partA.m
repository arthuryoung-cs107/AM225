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

prob2_part_a_15 = aysml_read('../dat_dir/prob2_Bruss_results_lambda-15');

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

delT_vec = zeros(size(precision_vec));

delT_vec(1)= 1e2;
delT_vec(2)= 5e2;
delT_vec(3)= 1e3;
delT_vec(4)= 5e3;
delT_vec(5)= 1e4;
delT_vec(6)= 5e4;
delT_vec(7)= 1e5;
delT_vec(8) = 5e5;
delT_vec(9) = 1e6;
delT_vec(10) = 5e6;
delT_vec(11) = 1e7;

delT_vec = (20.0)./( delT_vec );

figure(fig1.Number)
plot(prob6_parta_Brus_1e7(:, 1), prob6_parta_Brus_1e7(:, 2), '- ', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', 'y1: delT = 1e-1')
plot(prob6_parta_Brus_1e7(:, 1), prob6_parta_Brus_1e7(:, 3), '- ', 'Color', blue5, 'LineWidth', 1.5, 'DisplayName', 'y2: lamba = 1e-1')

figure(fig2.Number)
plot(delT_vec, precision_vec, '- ', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', 'convergence')
