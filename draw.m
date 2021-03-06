% drawing snr vs BLER(BER)

maxSNR = 8;
minSNR = -2;
snr = minSNR:0.05:maxSNR; 


% CBS = 528, Z = 64
%{
% -1.5~-0.2
CR_033  = [0.338 0.297 0.097 0.0716 0.0428 0.0225 0.00996 0.00832 0.00265 0.001385 0.000425 0.0001544 7.024e-5 4.45e-5];
% -1.1 -0.8 -0.5 -0.2
CR_033_mtk = [6.50E-02 1.00E-02 1.30E-03 1.50E-04];


% -0.8~1
CR_040  = [ 0.65 0.526 0.449 0.312  0.212 0.22 0.1 0.09 0.0547 0.0487 0.0172 0.00536 0.00312 0.001 0.000544 0.000466 0.000199 0.000124 7.69e-5];
% -0.05 0.25 0.55 0.85
CR_040_mtk = [5.50E-01 1.00E-02 1.30E-03 1.10E-04];


% 1.2~2.3
CR_050  = [0.16 0.068 0.0455 0.0185 0.011 0.00544 0.0021 0.0015 0.000478 0.000348 6.46e-5 6.04e-5];
% 1.35 1.65 1.95 2.25
CR_050_mtk = [5.50E-02 9.00E-03 1.20E-03 2.70E-04];



% 3.3~4.4
CR_067  = [0.146 0.102 0.0458 0.034 0.022 0.01 0.00432 0.00186 0.0013 0.000452 0.000302 0.000186];
% 3.45 3.7 4.05 4.4
CR_067_mtk = [5.00E-02 1.00E-02 1.50E-03 1.80E-04];



% 4.5~5.4
CR_075  = [0.107 0.0285 0.016 0.01 0.00577 0.00532 0.0014 0.0012 0.000521 0.000256];
%  4.55 4.85 5.15 5.45
CR_075_mtk = [5.00E-02 1.00E-02 1.20E-03 1.70E-04];



% 5.5~7
CR_083  = [0.15 0.157 0.077 0.0356 0.0228 0.0168 0.0112 0.005339 0.00536 0.002 0.00087 0.0006 0.00029 0.00013 9.47e-5 6.77e-5 ];
%  5.7 6 6.35 6.6 6.9
CR_083_mtk = [6.50E-02 1.80E-02 2.50E-03 6.20E-04 1.00E-04];



% 6.5~8
CR_089  = [0.2207 0.132 0.0644 0.0595 0.0483 0.0195 0.0153 0.0116 0.0066 0.00418 0.0026 0.00174 0.00104 0.000726 0.00034 0.0002126];
% 6.71 7.03 7.31 7.6
CR_089_mtk = [4.70E-02 1.20E-02 1.50E-03 2.20E-04];



% draw BLER
figure
semilogy(snr(11:2:37),CR_033,'k-diamond', snr([19,25,31,37]),CR_033_mtk,'r-diamond',...
         snr(25:2:61), CR_040,'k-diamond', snr([40,46,52,58]),CR_040_mtk,'r-diamond',...
         snr(65:2:87),CR_050,'k-diamond', snr([68,74,80,86]),CR_050_mtk,'r-diamond', ...
         snr(107:2:129), CR_067,'k-diamond', snr([110,115,122,129]),CR_067_mtk,'r-diamond',...
         snr(131:2:149),CR_075,'k-diamond', snr([132,138,144,150]),CR_075_mtk,'r-diamond',...
         snr(151:2:181), CR_083,'k-diamond', snr([155,161,168,173,179]),CR_083_mtk,'r-diamond', ...
         snr(171:2:201),CR_089,'k-diamond',snr([175,182,187,193]),CR_089_mtk,'r-diamond',  ...
         'LineWidth',1);
     

legend('code rate 0.33', 'code rate 0.33 mtk', ...
        'code rate 0.40', 'code rate 0.40 mtk', ...
        'code rate 0.50', 'code rate 0.50 mtk', ...
        'code rate 0.67', 'code rate 0.67 mtk', ...
        'code rate 0.75', 'code rate 0.75 mtk', ...
        'code rate 0.83', 'code rate 0.83 mtk', ...
        'code rate 0.89', 'code rate 0.89 mtk')

axis([minSNR, maxSNR+4, 10^-4, 1]);
hold on ; grid on;
xlabel('SNR (dB)');
ylabel('BLER');
title('BLER under different CR (CBS = 528)');
%}


% CBS = 1088, Z = 128
%{
% -1.4~-0.8
CR_033  = [0.064 0.0297 0.0148 0.00404 0.000808 0.00055 0.0000167];
% -1.3 -1.15 -1 -0.85 -0.7
CR_033_mtk = [3.10E-02 1.00E-02 1.95E-03 3.90E-04 1.00E-04];


% -0.5~0.3
CR_040  = [0.227 0.153 0.0278 0.00752 0.00433 0.00309 0.00124 0.000513 0.000197];
% -0.25 -0.1 0.1 0.25 0.4
CR_040_mtk = [3.40E-02 8.40E-03 1.60E-03 2.00E-04 1.00E-04];



% 1~1.8
CR_050  = [0.154 0.0445 0.0143 0.00702 0.002614 0.00165 0.000474 0.000125 9.749e-5];
% 1.2 1.35 1.5 1.6 1.75 1.85
CR_050_mtk = [2.30E-02 7.00E-03 1.60E-03 5.00E-04 2.00E-04 1.00E-04];



% 3~4.2
CR_067  = [0.35 0.117 0.0782 0.0353 0.0167 0.00892 0.00385 0.00175 0.000631 0.000567 0.000326 0.000268 0.000129];
% 3.3 3.45 3.6 3.75 3.9
CR_067_mtk = [2.40E-02 6.00E-03 1.40E-03 2.60E-04 1.00E-04];



% 4.4~5.2
CR_075  = [0.024 0.00798 0.00543 0.00266 0.00143 0.0009979 0.00053 0.000374 0.000214];
% 4.35 4.5 4.65 4.8 4.95
CR_075_mtk = [2.90E-02 8.40E-03 2.00E-03 4.90E-04 1.40E-04];



% 5.5~6.3
CR_083  = [0.053 0.022 0.00753 0.00401 0.00255 0.000815 0.000459 0.000224 0.000182];
%  5.5 5.6 5.8 5.95 6.1 6.3 6.35
CR_083_mtk = [6.00E-02 2.30E-02 7.80E-03 2.50E-03 6.90E-04 2.05E-04 1.00E-04];



% 6.5~7.5
CR_089  = [0.0888 0.0389 0.0345 0.0195 0.0205 0.00969 0.005 0.0022 0.000891 0.00071 0.00033];
% 6.55 6.7 6.85 7 7.15 7.3
CR_089_mtk = [2.80E-02 9.50E-03 3.10E-03 9.00E-04 3.10E-04 1.20E-04];


% draw BLER
figure
semilogy(snr(13:2:25),CR_033,'k-diamond', snr([15,18,21,24 27]),CR_033_mtk,'r-diamond', ...
         snr(31:2:47), CR_040,'k-diamond', snr([36,39,43,46,49]),CR_040_mtk,'r-diamond',...
         snr(61:2:77),CR_050,'k-diamond', snr([65,68,71,73,76,78]),CR_050_mtk,'r-diamond', ...
         snr(101:2:125), CR_067,'k-diamond', snr([107,110,113,116,119]),CR_067_mtk,'r-diamond',...
         snr(129:2:145),CR_075,'k-diamond', snr([128,131,134,137,140]),CR_075_mtk,'r-diamond',...
         snr(151:2:167), CR_083,'k-diamond', snr([151,153,157,160,163,167,168]),CR_083_mtk,'r-diamond',...
         snr(171:2:191),CR_089,'k-diamond',snr([172,175,178,181,184,187]),CR_089_mtk,'r-diamond',  ...
         'LineWidth',1);
legend('code rate 0.33', 'code rate 0.33 mtk', ...
        'code rate 0.40', 'code rate 0.40 mtk', ...
        'code rate 0.50', 'code rate 0.50 mtk', ...
        'code rate 0.67', 'code rate 0.67 mtk', ...
        'code rate 0.75', 'code rate 0.75 mtk', ...
        'code rate 0.83', 'code rate 0.83 mtk', ...
        'code rate 0.89', 'code rate 0.89 mtk')
axis([minSNR, maxSNR+4, 10^-4, 1]);
hold on ; grid on;
xlabel('SNR (dB)');
ylabel('BLER');
title('BLER under different CR (CBS = 1088)');
%}



% CBS = 2624, Z = 256
%{
% -1.6~-1.2
CR_033  = [0.0782 0.0329 0.00573 0.000933 0.000135];
% -1.6 -1.45 -1.3 -1.2
CR_033_mtk = [9.40E-02 1.40E-02 1.40E-03 1.00E-04];


% -0.6~-0.2
CR_040  = [0.27 0.0743 0.023 0.0015 0.0001];
% -0.5 -0.4 -0.25 -0.15
CR_040_mtk = [8.40E-02 1.20E-02 1.10E-03 1.00E-04];



% 0.8~1.3
CR_050  = [0.2245 0.061 0.0275 0.0046652 0.000862 0.00014];
% 0.85 1 1.15 1.3
CR_050_mtk = [5.40E-02 6.70E-03 5.40E-04 1.00E-04];



% 2.9~3.4
CR_067  = [0.33 0.0961 0.03386 0.0135 0.00151 0.000373  ];
% 2.95 3.1 3.25 3.4
CR_067_mtk = [7.90E-02 1.20E-02 9.80E-04 1.00E-04];



% 4~4.6
CR_075  = [0.104 0.029 0.018 0.00276 0.0011 0.0006625 0.000434];
%  4.05 4.2 4.35 4.5
CR_075_mtk = [5.90E-02 1.00E-02 9.50E-04 1.00E-04];



% 5.2~5.8
CR_083  = [0.101 0.033 0.00679 0.0029 0.000528 0.000176 5.14e-5];
%  5.25 5.4 5.55 5.7
CR_083_mtk = [4.30E-02 7.80E-03 8.70E-04 1.10E-04];



% 6.2~6.9
CR_089  = [0.168 0.0633 0.0457 0.0136 0.00685 0.00186 0.000941 0.0005];
% 6.2 6.4 6.5 6.65 6.8
CR_089_mtk = [5.60E-02 1.10E-02 1.90E-03 2.70E-04 1.00E-04];


% draw BLER
figure
semilogy(snr(9:2:17),CR_033,'k-diamond', snr([9,12,15,17]),CR_033_mtk,'r-diamond', ...
         snr(29:2:37), CR_040,'k-diamond', snr([31,33,36,38]),CR_040_mtk,'r-diamond',...
         snr(57:2:67),CR_050,'k-diamond', snr([58,61,64,67]),CR_050_mtk,'r-diamond', ...
         snr(99:2:109), CR_067,'k-diamond', snr([100,103,106,109]),CR_067_mtk,'r-diamond',...
         snr(121:2:133),CR_075,'k-diamond', snr([122,125,128,131]),CR_075_mtk,'r-diamond',...
         snr(145:2:157), CR_083,'k-diamond', snr([146,149,152,155]),CR_083_mtk,'r-diamond', ...
         snr(165:2:179),CR_089,'k-diamond',snr([165,169,171,174,177]),CR_089_mtk,'r-diamond',  ...
         'LineWidth',1);
legend('code rate 0.33', 'code rate 0.33 mtk', ...
        'code rate 0.40', 'code rate 0.40 mtk', ...
        'code rate 0.50', 'code rate 0.50 mtk', ...
        'code rate 0.67', 'code rate 0.67 mtk', ...
        'code rate 0.75', 'code rate 0.75 mtk', ...
        'code rate 0.83', 'code rate 0.83 mtk', ...
        'code rate 0.89', 'code rate 0.89 mtk')
axis([minSNR, maxSNR+4, 10^-4, 1]);
hold on ; grid on;
xlabel('SNR (dB)');
ylabel('BLER');
title('BLER under different CR (CBS = 2624)');
%}



% CBS = 4160, Z = 384
%{
% -1.7~-1.4
CR_033  = [0.09 0.0122 0.00147 0.000114];
% -1.65, -1.55, -1.45, -1.35, -1.3
CR_033_mtk = [8e-02 1.4e-2 1.5e-3 1.55e-4 1e-4];


% -0.7~-0.3
CR_040  = [0.385 0.103 0.0021 0.00247 0.000233];
% -0.6, -0.5, -0.4, -0.3
CR_040_mtk = [7.3e-2 1.6e-2 1.9e-3 1.4e-4];



% 0.7~1.1
CR_050  = [0.396 0.064 0.022 0.0026 0.000182];
% 0.85, 0.95, 1.05, 1.1
CR_050_mtk = [4.3e-2 5e-3 5.9e-4 1e-4];



% 2.9~3.3
CR_067  = [0.126 0.056 0.006 0.000645 6.28e-05 ];
% 2.95, 3, 3.1, 3.2, 3.25
CR_067_mtk = [6.5e-2 1.4e-2 2.05e-3 1.9e-4 1e-4];



% 3.9~4.4
CR_075  = [0.261 0.0695 0.011 0.0026 0.000575 0.0003079];
%  4, 4.1, 4.2,4.25, 4.3
CR_075_mtk = [5.4e-2 1.3e-2 1.7e-2 1.94e-4 1e-4];



% 5~5.5
CR_083  = [ 0.33 0.076 0.028 0.006 0.00156 0.000198 ];
%  5.15, 5.25, 5.35, 5.45, 5.55
CR_083_mtk = [4.8e-2 1.2e-2 1.8e-3 2.6e-4 1e-4];



% 6~6.8
CR_089  = [0.5 0.21 0.084 0.0304 0.0115 0.0026 7e-4 4e-4 9.9e-5];
% 6.2, 6.3, 6.4, 6.5, 6.6
CR_089_mtk = [3.3e-2 9e-3 2e-3 3.6e-4 1.1e-4];


% draw BLER
figure
semilogy(snr(7:2:13),CR_033,'k-diamond', snr([8,10,12,14,15]),CR_033_mtk,'r-diamond', ...
         snr(27:2:35), CR_040,'k-diamond', snr([29,31,33,35]),CR_040_mtk,'r-diamond',...
         snr(55:2:63),CR_050,'k-diamond', snr([58,60,62,63]),CR_050_mtk,'r-diamond', ...
         snr(99:2:107), CR_067,'k-diamond', snr([100,101,103,105,106]),CR_067_mtk,'r-diamond',...
         snr(119:2:129),CR_075,'k-diamond', snr([121,123,125,126,127]),CR_075_mtk,'r-diamond',...
         snr(141:2:151), CR_083,'k-diamond', snr([144,146,148,150,152]),CR_083_mtk,'r-diamond', ...
         snr(161:2:177),CR_089,'k-diamond',snr([165,167,169,171,173]),CR_089_mtk,'r-diamond',  ...
         'LineWidth',1);
legend('code rate 0.33', 'code rate 0.33 mtk', ...
        'code rate 0.40', 'code rate 0.40 mtk', ...
        'code rate 0.50', 'code rate 0.50 mtk', ...
        'code rate 0.67', 'code rate 0.67 mtk', ...
        'code rate 0.75', 'code rate 0.75 mtk', ...
        'code rate 0.83', 'code rate 0.83 mtk', ...
        'code rate 0.89', 'code rate 0.89 mtk')
axis([minSNR, maxSNR+2, 10^-4, 1]);
hold on ; grid on;
xlabel('SNR (dB)');
ylabel('BLER');
title('BLER under different CR (CBS = 4160)');
%}



% CBS = 7168, Z = 512
%{
% -1.8~-1.5
CR_033  = [0.2186 0.0283 0.00164 8.11e-05];
% -1.7, -1.65, -1.5, -1.45
CR_033_mtk = [1.8e-01 3e-02 2.1e-3 1e-4];


% -0.7~-0.4
CR_040  = [0.2132 0.0427 0.00598 0.000193];
% -0.65, -0.55, -0.45
CR_040_mtk = [9.9e-02 1.06e-02 4.9e-4];



% 0.6~1.0
CR_050  = [0.6176 0.1329 0.00832 0.000897 4.108e-5];
% 0.75, 0.85, 0.95, 1
CR_050_mtk = [1.2e-01 1.7e-02 1.5e-3 1e-4];



% 2.8~3.1
CR_067  = [0.326 0.095 0.0101 0.000564 ];
% 2.9, 3, 3.05, 3.1
CR_067_mtk = [6.3e-02 6.7e-03 3.3e-4 1e-4];



% 3.8~4.2
CR_075  = [0.456 0.089 0.0103 0.00103 5.25e-05];
% 3.9 , 4, 4.1, 4.2
CR_075_mtk = [7e-02 1.1e-02 6.4e-4 1e-4];



% 5~5.4
CR_083  = [ 0.101 0.023 0.00298 0.000339 0.00015 ];
% 5.05, 5.15, 5.25, 5.35
CR_083_mtk = [7e-02 1.2e-02 1.1e-3 1e-4];



% 6~6.4
CR_089  = [0.3125 0.0854 0.0162 0.00329 0.00023];
% 6.05, 6.15, 6.25, 6.35, 6.4
CR_089_mtk = [7.01e-02 1.6e-02 1.7e-3 1.7e-4 1e-4];



% draw BLER
figure
semilogy(snr(5:2:11),CR_033,'k-diamond', snr([7,8,11,12]),CR_033_mtk,'r-diamond', ...
         snr(27:2:33), CR_040,'k-diamond', snr(28:2:32),CR_040_mtk,'r-diamond',...
         snr(53:2:61),CR_050,'k-diamond', snr([58,60,62,63]),CR_050_mtk,'r-diamond', ...
         snr(97:2:103), CR_067,'k-diamond', snr([99,101,102,103]),CR_067_mtk,'r-diamond',...
         snr(117:2:125),CR_075,'k-diamond', snr([119,121,123,125]),CR_075_mtk,'r-diamond',...
         snr(141:2:149), CR_083,'k-diamond', snr([142,144,146,148]),CR_083_mtk,'r-diamond', ...
         snr(161:2:169),CR_089,'k-diamond',snr([162,164,166,168,169]),CR_089_mtk,'r-diamond',  ...
         'LineWidth',1);
legend('code rate 0.33', 'code rate 0.33 mtk', ...
        'code rate 0.40', 'code rate 0.40 mtk', ...
        'code rate 0.50', 'code rate 0.50 mtk', ...
        'code rate 0.67', 'code rate 0.67 mtk', ...
        'code rate 0.75', 'code rate 0.75 mtk', ...
        'code rate 0.83', 'code rate 0.83 mtk', ...
        'code rate 0.89', 'code rate 0.89 mtk')
axis([minSNR, maxSNR+2, 10^-4, 1]);
hold on ; grid on;
xlabel('SNR (dB)');
ylabel('BLER');
title('BLER under different CR (CBS = 7168)');
%}



% CBS = 7808, Z = 512
%{
% -1.8~-1.5
CR_033  = [0.3499 0.0589 0.0031 16e-05];
% -1.75, -1.65, -1.55
CR_033_mtk = [1.2e-01 1.3e-02 4.4e-4];


% -0.7~-0.4
CR_040  = [0.1973 0.0227 0.0019 4e-05];
% -0.65, -0.55, -0.45
CR_040_mtk = [7.5e-02 7.5e-03 3.1e-4];



% 0.6~0.9
CR_050  = [0.4394 0.1001 0.0086 2.7548e-04];
% 0.8, 0.9, 0.95
CR_050_mtk = [4e-02 6e-04 1e-4];



% 2.8~3.1
CR_067  = [0.1193 0.0171 0.0015 7e-05 ];
% 2.85, 2.95, 3.05, 3.1
CR_067_mtk = [6.5e-02 8e-03 3.6e-4 1e-4];



% 3.8~4.1
CR_075  = [0.3214 0.0669 0.0085 5.8634e-04];
% 3.9 , 4, 4.05, 4.15
CR_075_mtk = [1.05e-01 1.04e-02 1.01e-3 1e-4];



% 5~5.3
CR_083  = [ 0.1711 0.0369 0.0055 5.1243e-04 ];
% 5.05, 5.15, 5.25, 5.35
CR_083_mtk = [9.4e-02 1.4e-02 1.4e-3 1.3e-4];



% 6~6.4
CR_089  = [0.1672 0.0349 0.0059 4.6512e-04 9e-05];
% 6.05, 6.15, 6.25, 6.35, 6.45
CR_089_mtk = [8.1e-02 1.8e-02 2e-3 2.1e-4 1e-4];



% draw BLER
figure
semilogy(snr(5:2:11),CR_033,'k-diamond', snr(6:2:10),CR_033_mtk,'r-diamond', ...
         snr(26:2:32), CR_040,'k-diamond', snr(27:2:31),CR_040_mtk,'r-diamond',...
         snr(53:2:59),CR_050,'k-diamond', snr([57,59,60]),CR_050_mtk,'r-diamond', ...
         snr(89:2:95), CR_067,'k-diamond', snr(90:2:96),CR_067_mtk,'r-diamond',...
         snr(109:2:115),CR_075,'k-diamond', snr([111,113,114,116]),CR_075_mtk,'r-diamond',...
         snr(132:2:138), CR_083,'k-diamond', snr(133:2:139),CR_083_mtk,'r-diamond', ...
         snr(153:2:161),CR_089,'k-diamond',snr(154:2:162),CR_089_mtk,'r-diamond',  ...
         'LineWidth',1);
legend('code rate 0.33', 'code rate 0.33 mtk', ...
        'code rate 0.40', 'code rate 0.40 mtk', ...
        'code rate 0.50', 'code rate 0.50 mtk', ...
        'code rate 0.67', 'code rate 0.67 mtk', ...
        'code rate 0.75', 'code rate 0.75 mtk', ...
        'code rate 0.83', 'code rate 0.83 mtk', ...
        'code rate 0.89', 'code rate 0.89 mtk')
axis([minSNR, maxSNR+2, 10^-4, 1]);
hold on ; grid on;
xlabel('SNR (dB)');
ylabel('BLER');
title('BLER under different CR (CBS = 7808)');
%}




