function errorStats = powerDelayProf(channel)
% txData = randi([0 1],4096*60418,1);
load('txData.mat');
bpskModulator = comm.BPSKModulator;
bpskModulator.PhaseOffset = 0;
modData = bpskModulator(txData);
modDataFreq = fft(modData,length(channel));
recDataFreq = channel'.*modDataFreq;
recData = ifft(recDataFreq,4096*60418);
bpskDemodulator = comm.BPSKDemodulator;
errorRate = comm.ErrorRate;
rxData = bpskDemodulator(recData);       % Demodulate
pckLoss=0;
for i = 1:60418
    errorStats = 0;
    errorStats = biterr(txData((1+(i-1)*(4096)):4096*i),rxData((1+(i-1)*(4096)):4096*i)); % Collect error stats
    if errorStats > 2112
        pckLoss = pckLoss+1;
    end
end
%bitErrRate = biterr(txData,rxData);
end
