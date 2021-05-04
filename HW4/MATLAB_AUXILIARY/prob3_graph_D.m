fig1 = figure('Name', 'alt cubic solution', 'Renderer', 'painters', 'Position', fig_pos(1, :));
xlabel('x')
ylabel('u_{fea}')
hold on

fig2 = figure('Name', 'analytical solution', 'Renderer', 'painters', 'Position', fig_pos(6, :));
xlabel('x')
ylabel('u_{sol}')
hold on

fig3 = figure('Name', 'C2 residual', 'Renderer', 'painters', 'Position', fig_pos(2, :));

N_test=1000;

N_min = 10;
delta = 5;
N_max = 1000;

u_sol_fun = @(x) exp(1-x).*sin(5*pi*x);

figure(fig2.Number)
plot(1:0.001:2 , u_sol_fun(1:0.001:2), ' -', 'Color', blue5, 'LineWidth', 1.5, 'DisplayName', 'analytical solution')

figure(fig1.Number)
for N=N_min:delta:N_max
  clf(fig1)
  clf(fig3)

  prob3_sol = aysml_read(['../dat_dir/prob3_altcube_C1_N', num2str(N),'_Ntest', num2str(N_test)]);
  u_sol = u_sol_fun(prob3_sol(:, 1)) ;

  figure(fig1.Number)
  xlabel('x')
  ylabel('u_{fea}')
  hold on
  plot(prob3_sol(:, 1), prob3_sol(:, 2), ' o', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', 'FEA solution')
  plot(prob3_sol(:, 1), u_sol, ' o', 'Color', blue5, 'LineWidth', 1.5, 'DisplayName', 'analytical solution')
  legend('Show')

  figure(fig3.Number);
  xlabel('x')
  ylabel('residual')
  hold on
  plot(prob3_sol(:, 1), u_sol-prob3_sol(:, 2), ' o', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', 'C2 FEA solution')

  N
  pause(1);

end
