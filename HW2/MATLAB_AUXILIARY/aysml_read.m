function mat_return = aysml_read(name)
  dims = dlmread([name '.aysml']);
  m = dims(1);
  n = dims(2);

  mat_return = (fread(fopen([name '.aydat']), [n, m], 'float64=>float64'))';

end
