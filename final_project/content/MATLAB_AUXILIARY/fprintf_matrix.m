function fprintf_matrix(mat, name)
  file_id = fopen([name, '.aydat'], 'w+');

  for i=1:size(mat, 1)
    for j=1:size(mat, 2)
      fwrite(file_id, mat(i, j), 'double');
    end
  end
  fclose(file_id);
  file_id2 = fopen([name, '.aysml'], 'w+');
  fprintf(file_id2, '%d %d', size(mat, 1), size(mat, 2));
  fclose(file_id2);
end
