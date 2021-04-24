fig1 = figure('Name', 'l2 error vs. N', 'Renderer', 'painters', 'Position', fig_pos(1, :));
xlabel('')
ylabel('')
hold on

x_vec = 0:0.01:1;
phi_vec = zeros(size(x_vec));

for i=1:length(x_vec)
  phi_vec(i) = phi_cubic(x_vec(i));
end

figure(fig1.Number)
plot( x_vec, phi_vec, ' o', 'Color', red4, 'LineWidth', 1.5, 'DisplayName', '')
legend('Show')
