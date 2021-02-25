snr = 0:5:30;
bw = 0.08E9;
seProb = 10E-5;
kn = 4;

[channel, rad_res_TE, rad_res_TM, rad_reac_TE, rad_reac_TM, ant_res_TE, ant_res_TM, ant_reac_TE, ant_reac_TM, ant_imp_TE, ant_imp_TM, TE_mode_imp, TM_mode_imp] = chImpRespCyl_FreqBand();
radRes = [rad_res_TE rad_res_TM];
radReac = [rad_reac_TE rad_reac_TM];
radImp = [TE_mode_imp TM_mode_imp];
antRes = [ant_res_TE ant_res_TM];
antReac = [ant_reac_TE ant_reac_TM];
antImp = [ant_imp_TE ant_imp_TM];

for i = 1:length(snr)
    rate(i) = bw*log2(1 + ((3*snr(i)*((sum(channel))^2))/((erfcinv(seProb/kn))^2)));
end
figure
plot(snr,rate);
