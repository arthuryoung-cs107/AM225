fig1 = figure('Name', 'matrix mult times, generic method', 'Renderer', 'painters', 'Position', fig_pos(1, :));
xlabel('n')
ylabel('time (s)')
hold on

fig2 = figure('Name', 'matrix mult times, Strassen recursion', 'Renderer', 'painters', 'Position', fig_pos(2, :));
xlabel('n')
ylabel('time (s)')
hold on

fig3 = figure('Name', 'matrix mult times, BLAS', 'Renderer', 'painters', 'Position', fig_pos(3, :));
xlabel('n')
ylabel('time (s)')
hold on


timesA = dlmread('../dat_dir/prob1_partA_times8.dat');
timesB = dlmread('../dat_dir/prob1_partB_times8.dat');
timesC = dlmread('../dat_dir/prob1_partC_times8.dat');

fA = fit(timesA(:, 1),timesA(:, 2),'b*x^m');
fB = fit(timesB(:, 1),timesB(:, 2),'b*x^m');
fC = fit(timesC(:, 1),timesC(:, 2),'b*x^m');

figure(fig1.Number)
plot(timesA(:, 1), timesA(:, 2), '- s', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', 'generic')
fplot(@(n) fA.b.*(n).^(fA.m), [min(timesA(:, 1)), max(timesA(:, 1)) ], '- k', 'Linewidth', 2,'DisplayName', 'power law fit')
legend('Show', 'Location', 'NorthWest')
set(gca, 'XTickLabel',[])                      %# suppress current x-labels
xt = get(gca, 'XTick');
yl = get(gca, 'YLim');
str = cellstr( num2str(xt(:),'2^{%d}') );      %# format x-ticks as 2^{xx}
hTxt = text(xt, yl(ones(size(xt))), str, ...   %# create text at same locations
    'Interpreter','tex', ...                   %# specify tex interpreter
    'VerticalAlignment','top', ...             %# v-align to be underneath
    'HorizontalAlignment','center');           %# h-aligh to be centered

figure(fig2.Number)
plot(timesB(:, 1), timesB(:, 2), '- s', 'Color', blue5, 'LineWidth', 1.5, 'DisplayName', 'Strassen')
fplot(@(n) fB.b.*(n).^(fB.m), [min(timesA(:, 1)), max(timesA(:, 1)) ], '- k', 'Linewidth', 2,'DisplayName', 'power law fit')
legend('Show', 'Location', 'NorthWest')
set(gca, 'XTickLabel',[])                      %# suppress current x-labels
xt = get(gca, 'XTick');
yl = get(gca, 'YLim');
str = cellstr( num2str(xt(:),'2^{%d}') );      %# format x-ticks as 2^{xx}
hTxt = text(xt, yl(ones(size(xt))), str, ...   %# create text at same locations
    'Interpreter','tex', ...                   %# specify tex interpreter
    'VerticalAlignment','top', ...             %# v-align to be underneath
    'HorizontalAlignment','center');           %# h-aligh to be centered

figure(fig3.Number)
plot(timesC(:, 1), timesC(:, 2), '- s', 'Color', green4, 'LineWidth', 1.5, 'DisplayName', 'BLAS')
fplot(@(n) fC.b.*(n).^(fC.m), [min(timesA(:, 1)), max(timesA(:, 1)) ], '- k', 'Linewidth', 2,'DisplayName', 'power law fit')
legend('Show', 'Location', 'NorthWest')
set(gca, 'XTickLabel',[])                      %# suppress current x-labels
xt = get(gca, 'XTick');
yl = get(gca, 'YLim');
str = cellstr( num2str(xt(:),'2^{%d}') );      %# format x-ticks as 2^{xx}
hTxt = text(xt, yl(ones(size(xt))), str, ...   %# create text at same locations
    'Interpreter','tex', ...                   %# specify tex interpreter
    'VerticalAlignment','top', ...             %# v-align to be underneath
    'HorizontalAlignment','center');           %# h-aligh to be centered
