prefix = '../aydat_dir_sim2/';
location = '../../../../sim2_results/simh_3_620K_10K/';

for i=0:500
  rydat2aydat([location, 'tem.', num2str(i)], [prefix, 'tem', num2str(i)]);
end
