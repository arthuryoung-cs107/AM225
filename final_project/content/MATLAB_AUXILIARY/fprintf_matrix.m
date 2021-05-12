function fprintf_matrix(mat, name)
  file_id = fopen([name, '.aydat'], 'w+');
  mat_trans = mat';
  fwrite(file_id, mat_trans(:), 'double');
  fclose(file_id);
  file_id2 = fopen([name, '.aysml'], 'w+');
  fprintf(file_id2, '%d %d', size(mat, 1), size(mat, 2));
  fclose(file_id2);
end
