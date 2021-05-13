function rydat2aydat_corrupt(loc, prefix)
  s = rydat_read(loc);
  rand_mat = (unifrnd(0, 1, size(s.mat))) < 0.1; % percent sparsity
  bin_rand = 2*(randi([0, 1], size(s.mat))) - 1;
  x = s.mat + 10*std(s.mat(:))*(rand_mat.*bin_rand);
  fprintf_matrix(x, prefix);
end
