function struct_return = rydat_read(name)
  id = fopen(name);
  mat_read = (fread(id, 'float32'))';
  fclose(id);
  x_length = mat_read(1);
  y_length = (length(mat_read)- (x_length +1))/(x_length + 1);

  x_vec = mat_read(2:(x_length+1));
  y_vec = zeros(1, y_length);
  mat_out = zeros(y_length, x_length);

  index = 1 + x_length;
  for i=1:y_length
    index = index + 1;
    y_vec(i) = mat_read(index);
    for j=1:x_length
      index = index + 1;
      mat_out(i, j) = mat_read(index);
    end
  end

  struct_return.x = x_vec;
  struct_return.y = y_vec;
  struct_return.mat = mat_out;
end
