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

dims = dlmread('../dat_dir/prob2_specs_threads8.dat');
m = dims(1);
n = dims(2);

fig1 = figure('Name', 'Jet space, view 1', 'Renderer', 'painters', 'Position', fig_pos(1, :));
ylabel('efficiency')
% xlabel('n')
box on
hold on

% un-noised G matrix
thread8_results = (fread(fopen('../dat_dir/prob2_results_threads8.dat'), [n, m], 'float64=>float64'))';
thread1_results = (fread(fopen('../dat_dir/prob2_results_threads1.dat'), [n, m], 'float64=>float64'))';
thread2_results = (fread(fopen('../dat_dir/prob2_results_threads2.dat'), [n, m], 'float64=>float64'))';
thread3_results = (fread(fopen('../dat_dir/prob2_results_threads3.dat'), [n, m], 'float64=>float64'))';
thread4_results = (fread(fopen('../dat_dir/prob2_results_threads4.dat'), [n, m], 'float64=>float64'))';

thread_all_results = [thread1_results; thread2_results; thread3_results; thread4_results];
n_all = thread_all_results(:, 2);
n_all = reshape(n_all,[m, 4]);
time_all = thread_all_results(:, 3);
time_all = reshape(time_all,[m, 4]);

efficiency = zeros(7, 4);

for i = 1:m
  for j = 1:4
    efficiency(i, j) = time_all(i, 1)/( j*time_all(i, j));
  end
end

figure(fig1.Number)
plot(log2(n_all(:, 2)), efficiency(:, 2), '- s', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', '2 threads')
plot(log2(n_all(:, 3)), efficiency(:, 3), '- o', 'Color', blue5, 'LineWidth', 1.5, 'DisplayName', '3 threads')
plot(log2(n_all(:, 4)), efficiency(:, 4), '- ^', 'Color', green4, 'LineWidth', 1.5, 'DisplayName', '4 threads')
legend('Show', 'Location', 'SouthEast')
set(gca, 'XTickLabel',[])                      %# suppress current x-labels
xt = get(gca, 'XTick');
yl = get(gca, 'YLim');
str = cellstr( num2str(xt(:),'2^{%d}') );      %# format x-ticks as 2^{xx}
hTxt = text(xt, yl(ones(size(xt))), str, ...   %# create text at same locations
    'Interpreter','tex', ...                   %# specify tex interpreter
    'VerticalAlignment','top', ...             %# v-align to be underneath
    'HorizontalAlignment','center');           %# h-aligh to be centered
