function out = phi_cubic(x_in, i) %% notice: i is index 1
  omega = [1, 2];
  N = 3;
  N_full = 3*N + 1;

  if (x_in==omega(2) && i==N_full)
    out = 1;
  else
    h = (omega(2)-omega(1))/(N_full-1);
    node_pos = (omega(1)):h:(omega(2));
    i_mod = mod((i-1), 3);

    if (i_mod == 0 ) %% edge node?
      if ( abs(x_in - node_pos(i))  < h ) %% in proximity of edge node
        if (x_in >= node_pos(i)) %% in primary interval
          out = 1;
          for k=1:1:3
            out = out*(x_in - node_pos(i + k))/(node_pos(i) - node_pos(i + k));
          end
        else
          out = 1;
          for k=-3:1:-1
            out = out*(x_in - node_pos(i + k))/(node_pos(i) - node_pos(i + k));
          end
        end
      else
        out = 0;
      end
    else %% is an interior node
      if ( (node_pos(i-i_mod) < x_in) && (node_pos(i-i_mod + 3) > x_in ) ) %% if in the relevant interval
        out = 1;
        for k= (i-i_mod):(i-1)
          out = out*(x_in - node_pos(k))/(node_pos(i) - node_pos(k));
        end
        for k=(i+1):(i-i_mod + 3)
          out = out*(x_in - node_pos(k))/(node_pos(i) - node_pos(k));
        end
      else
        out = 0;
      end
    end
end
