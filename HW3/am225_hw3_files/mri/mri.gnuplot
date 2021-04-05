set terminal png
set output 'mri_data.png'
unset xtics
unset ytics
set yrange [*:*] reverse
plot 'out' matrix with image 
