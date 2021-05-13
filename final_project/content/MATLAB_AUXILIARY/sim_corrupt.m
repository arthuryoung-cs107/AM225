prefix = '../aydat_dir_small_sim5/';
location = '../../../../sim5_small_results/simh_3_620K_10K/';

for i=0:500
  rydat2aydat_corrupt([location, 'tem.', num2str(i)], [prefix, 'tem_corrupt', num2str(i)]);
end
