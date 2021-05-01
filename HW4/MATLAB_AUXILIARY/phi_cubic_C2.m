function out = phi_cubic_C2(x)

  if (-2 <= x && x <= -1)
    out = 0.25*(x + 2)^3;
  elseif (-1 <= x && x <= 0)
    out = 0.25*( 1 + 3*(1 + x) + 3*(1 + x)^2 - 3*(1 + x)^3 );
  elseif (0 <= x && x <= 1)
    out = 0.25*( 1 + 3*(1 - x) + 3*(1 - x)^2 - 3*(1 - x)^3 );
  elseif (1 <= x && x <= 2)
    out = 0.25*(2-x)^3;
  end

end
