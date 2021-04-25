function out = phi_cubic_C1(x_in, i) %% notice: i is index 1
  omega = [1, 2];
  N = 3;
  N_full = N + 1;
  h = (omega(2) - omega(1))/(N);
  node_pos = (omega(1)):(h):(omega(2));

  x_eval = abs(x_in - node_pos(i));

  if (x_eval < h)
    out = 1 - 3*((x_eval/h)^2) + 2*((x_eval/h)^3);
  else
    out = 0;
end
