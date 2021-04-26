fig1 = figure('Name', 'alt cubic solution', 'Renderer', 'painters', 'Position', fig_pos(1, :));
xlabel('x')
ylabel('u_{fea}')
hold on

fig2 = figure('Name', 'analytical solution', 'Renderer', 'painters', 'Position', fig_pos(2, :));
xlabel('x')
ylabel('u_{sol}')
hold on

N_test=1000;

N_min = 10;
delta = 10;
N_max = 1000;

u_sol_fun = @(x) exp(1-x).*sin(5*pi*x);

figure(fig1.Number)
for N=N_min:delta:N_max
  clf(fig1)
  xlabel('x')
  ylabel('u_{fea}')
  hold on

  prob3_sol = aysml_read(['../dat_dir/prob3_altcube_N', num2str(N),'_Ntest', num2str(N_test)]);
  u_sol = u_sol_fun(prob3_sol(:, 1)) ;
  plot(prob3_sol(:, 1), prob3_sol(:, 2), ' o', 'Color', red5, 'LineWidth', 1.5, 'DisplayName', 'FEA solution')
  plot(prob3_sol(:, 1), u_sol, ' o', 'Color', blue5, 'LineWidth', 1.5, 'DisplayName', 'analytical solution')
  legend('Show')

  pause(1);

end
