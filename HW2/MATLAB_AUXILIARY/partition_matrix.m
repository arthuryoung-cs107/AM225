function mat_out = partition_matrix(r, bush)
  mat_out = zeros(1, bush);
  mat_out(1, 1) = r;
  width = 2;
  row_last = mat_out(size(mat_out, 1), :);

  while (width <= bush && row_last(1) > 1)
    row_new = row_last;
    row_new(1) = row_last(1) - 1;
    row_new(width) = 1;
    mat_out = [mat_out; row_new];
    mat_out = bush_recurse(mat_out, width);

    width = width + 1;
    row_last = row_new;
  end
end
