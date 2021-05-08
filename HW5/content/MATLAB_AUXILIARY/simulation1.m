fig1 = figure('Name', 'simulation1', 'Renderer', 'painters', 'Position', fig_pos(1, :));
box on
hold on



for i = 0:199
  red = aysml_read(['../dat_dir/sim_colourR_frame', num2str(i)]);
  green = aysml_read(['../dat_dir/sim_colourG_frame', num2str(i)]);
  blue = aysml_read(['../dat_dir/sim_colourB_frame', num2str(i)]);

  RGB = zeros(size(red, 1), size(red, 2), 3);
  RGB(:, :, 1) = red;
  RGB(:, :, 2) = green;
  RGB(:, :, 3) = blue;
  
  clf(fig1);
  figure(fig1.Number);
  imshow(RGB);

  pause(0.01);

end

red = aysml_read(['../dat_dir/sim_colourR_frame', num2str(199)]);
green = aysml_read(['../dat_dir/sim_colourG_frame', num2str(199)]);
blue = aysml_read(['../dat_dir/sim_colourB_frame', num2str(199)]);

RGB = zeros(size(red, 1), size(red, 2), 3);
RGB(:, :, 1) = red;
RGB(:, :, 2) = green;
RGB(:, :, 3) = blue;


figure(fig1.Number);
imshow(RGB);
