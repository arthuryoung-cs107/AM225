% X_check = rydat_read('../../../../sim1_results/simh_3_620K_20K/X.0');
% Y_check = rydat_read('../../../../sim1_results/simh_3_620K_20K/Y.0');
% tem_check = rydat_read('../../../../sim1_results/simh_3_620K_20K/tem.0');
% p_check = rydat_read('../../../../sim1_results/simh_3_620K_20K/p.0');
% dev_check = rydat_read('../../../../sim1_results/simh_3_620K_20K/dev.0');
%
prefix = '../aydat_dir_sim1/';
location = '../../../../sim1_results/simh_3_620K_20K/';

for i=0:600
  % rydat2aydat([location, 'X.', num2str(i)], [prefix, 'X', num2str(i)]);
  % rydat2aydat([location, 'Y.', num2str(i)], [prefix, 'Y', num2str(i)]);
  rydat2aydat([location, 'tem.', num2str(i)], [prefix, 'tem', num2str(i)]);
  rydat2aydat([location, 'p.', num2str(i)], [prefix, 'p', num2str(i)]);
  rydat2aydat([location, 'dev.', num2str(i)], [prefix, 'dev', num2str(i)]);
end
