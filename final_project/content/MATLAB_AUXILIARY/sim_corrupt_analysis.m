clear
close all

fig_pos = [10, 455, 350, 350;
    360, 455, 350, 350; 710, 455, 350, 350; 1060, 455, 350, 350; 10, 30, 350, 350; 360, 30, 350, 350; 710, 30, 350, 350; 1060, 30, 350, 350];

grey1 = [178/250, 186/250, 187/250];
grey2 = [131/250, 145/250, 146/250];
grey3 = [97/250, 106/250, 107/250];
grey4 = [66/250, 73/250, 73/250];
grey5 = [20/100 20/100 20/100];
grey = [grey1; grey2; grey3; grey4; grey5];

purple1 = [102/250, 0/250, 102/250];
purple2 = [153/250, 0/250, 153/250];
purple3 = [204/250, 0/250, 204/250];
purple4 = [250/250, 0/250, 250/250];
purple5 = [250/250, 50/250, 250/250];
purple = [purple1; purple2; purple3; purple4; purple5];

orange1 = [255/255 90/255 0];
orange2 = [255/255 123/255 0];
orange3 = [255/255 165/255 0];
orange4 = [255/255 208/255 0];
orange5 = [255/255 229/255 0];
orange = [orange1; orange2; orange3; orange4; orange5];

green1 = [88/250, 214/250, 141/250];
green2 = [40/250, 180/250, 99/250];
green3 = [34/250, 153/250, 84/250];
green4 = [25/250, 111/250, 61/250];
green5 = [0, 1, 0];
green = [green1; green2; green3; green4; green5];

blue1 = [120/250, 150/250, 250/250];
blue2 = [52/250, 152/250, 219/250];
blue3 = [39/250, 97/250, 141/250];
blue4 = [10/250, 50/250, 150/250];
blue5 = [0, 0, 1];
blue = [blue1; blue2; blue3; blue4; blue5];

red1 = [236/250, 112/250, 99/250];
red2 = [192/250, 57/250, 43/250];
red3 = [146/250, 43/250, 33/250];
red4 = [100/250, 30/250 , 22/250];
red5 = [1, 0 , 0];
red = [red1; red2; red3; red4; red5];

pink = [255 0 104 ; 243 0 112; 230 0 119 ; 216  0 125; 200 0 131; 183 0 136; 165 0 140; 146 0 143; 125 7 145; 102 20 146]*(255)^(-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run sim3_corrupt_analysis.m
% run sim4_corrupt_analysis.m
% run sim5_corrupt_analysis.m

% run simALL_corrupt_analysis.m
% run simALL_PC_analysis.m
% run simALL_frame_by_frame.m
run simALL_transition.m
