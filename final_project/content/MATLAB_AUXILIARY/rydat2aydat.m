function rydat2aydat(loc, prefix)
  s = rydat_read(loc);
  fprintf_matrix(s.mat, prefix);
end
