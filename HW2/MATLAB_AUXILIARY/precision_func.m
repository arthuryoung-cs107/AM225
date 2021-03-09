function precision = precision_func(input, reference)
  dy0 = reference(size(reference, 1), 2) - input(size(input, 1), 2);
  dy1 = reference(size(reference, 1), 3) - input(size(input, 1), 3);
  precision = sqrt(dy0*dy0 + dy1*dy1);
end
