prefix = '../aydat_dir_small_sim4/';
location = '../../../../sim4_small_results/simh_3_620K_20K/';

for i=0:500
  rydat2aydat([location, 'tem.', num2str(i)], [prefix, 'tem', num2str(i)]);
end
