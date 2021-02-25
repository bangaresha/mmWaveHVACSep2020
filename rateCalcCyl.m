function rate = rateCalcCyl(bw,snr,channel,seProb,kn)
for i = 1:length(snr)
    rate(i) = bw*log2(1 + ((3*snr(i)*((sum(channel))^2))/((erfcinv(seProb/kn))^2)));
end
figure
plot(snr,rate);
end
