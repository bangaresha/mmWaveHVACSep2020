function rate = rateCalcCyl16QAM(bw,snr,channel,seProb,kn,fgap)

ko = ((erfcinv(seProb/kn))^2);
% for j = 1:length(snr)
%     for i = 1:length(channel)
% %         funclam = funclam + ((lamb*fgap/log(2)) - ...
% %             fgap/(bw*(3*snr(j)./((erfcinv(seProb/kn))^2))*(channel(i).^2)));
%           funclam = funclam + (fgap*ko)/(bw*3*(snr(j))*(channel(i).^2));
%     end
%     lambda(j,i) = (1 + funclam)*log(2)/fgap;
% end

for i = 1:length(snr)
    rateTemp=0;
    for j=1:length(channel)
        if channel(j) == 0
            channel(j) = -0.0001;
        end
        deno = (length(channel))*ko*log(2);
        nume = 3*bw*(snr(i))*(abs(channel(j)).^2);
        frac = nume/deno;
        rateTemp = rateTemp + fgap*log2(frac);
    end
    rate(i) = rateTemp/24;
end

figure
plot(snr,rate);
end

