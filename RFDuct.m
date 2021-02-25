distance = [1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5];
RSSI = [-36, -37, -38, -39, -39, -40, -41, -43];
SNR = [13, 12, 9, 8, 8, 7, 6, 5];
EVM = [-10, -9, -8, -8, -7, -7, -6, -5];
BitRate = [640.3, 638.9, 638, 637, 636, 635, 634, 632.8];

figure
yyaxis left
plot(distance, RSSI);
yyaxis right
plot(distance, BitRate);
title('RSSI and Bit Rate versus distance')
yticks('auto')
xticks('auto')

% xticks([0.81, 1.06, 1.37, 1.72, 2.03, 2.31, 2.61, 2.89, 3.14, 3.47, 3.75, 4.03, 4.29, 4.57, 4.85])
% xlim([0.81 4.85])

% figure
% plot(distance, SNR);
% title('SNR versus distance')
% yticks('auto')
% xticks('auto')

figure
yyaxis left
plot(distance, EVM);
yyaxis right
plot(distance, SNR);
%title('EVM and SNR versus distance')
yticks('auto')
xticks('auto')



% figure
% subplot(2,2,1)
% plot(distance, RSSI);
% title('RSSI versus distance')
% subplot(2,2,2)
% plot(distance, SNR);
% title('SNR versus distance')
% subplot(2,2,3)
% plot(distance, EVM);
% title('EVM versus distance')
% subplot(2,2,4)
% plot(distance, BitRate);
% title('Bit Rate versus distance')

%Coefficients (with 95% confidence bounds):
        a1 =      0.0823;%  (0.07829, 0.08632)
       b1 =   3.816e-10;%  (fixed at bound)
       c1 =       47.02;%  (46.96, 47.08)
       a2 =     0.02753;%  (0.02321, 0.03185)
       b2 =   5.175e-09;%  (fixed at bound)
       c2 =      -67.59;%  (-67.75, -67.43)
       freq = 59E9:4E6:60E9;
       omega = 2*pi*freq;
       f=[];
       for x = omega
           %f =  [f a1*sin(b1*x+c1)];
           f =  [f a1*sin(b1*x+c1) + a2*sin(b2*x+c2)];
       end
       h = ifft(f);
       figure
       plot(abs(f));
       figure
       plot(abs(h));
       %188.46 -0.0121   1850.82 0.7751 0.0213 0.00919