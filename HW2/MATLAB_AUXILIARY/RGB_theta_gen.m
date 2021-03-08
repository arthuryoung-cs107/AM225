function RGB_out = RGB_theta_gen(theta_in)
  f_theta = @(theta) 0.45*(1 + cos(theta));
  RGB_out = zeros(3, length(theta_in));

  RGB_out(1, :) = f_theta(theta_in);
  RGB_out(2, :) = f_theta(theta_in - 2*pi/3);
  RGB_out(3, :) = f_theta(theta_in + 2*pi/3);
  RGB_out = RGB_out';
end
