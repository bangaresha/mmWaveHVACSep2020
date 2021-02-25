function rate = rateCalcCylBpsk(bw,snr,channel,seProb,fgap)

for i = 1:length(snr)
    rateTemp=0;
    pb = 0;
    for j=1:length(channel)
        pb = 0.5*erfc(sqrt(bw*(10^((bw*snr(i))/10))));
        deno = (length(channel))*((erfcinv(seProb))^2)*log(2);
        nume = bw*(snr(i))*(channel(j).^2);
        frac = nume/deno;
        rateTemp = rateTemp + fgap*log2(frac);
    end
    rate(i) = rateTemp;
end

figure
plot(snr,rate);
end

