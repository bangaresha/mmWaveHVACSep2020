distance = [4; 8];
angle = [0, 45, 90, 180];
RSSI4 = [-37, -40, -44, -42];
RSSI8 = [-38, -44, -48, -45];
BitRate4 = [544, 542, 538, 540];
BitRate8 = [542, 539, 535, 540];

figure
subplot(2,1,1)
plot(angle, RSSI4);
subplot(2,1,2)
subplot(angle, BitRate4);
yticks('auto')
xticks('auto')
figure
subplot(2,1,1)
plot(angle, RSSI8);
subplot(2,1,2)
plot(angle, BitRate8);
yticks('auto')
xticks('auto')





