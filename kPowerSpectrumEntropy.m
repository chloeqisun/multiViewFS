function psdE = kPowerSpectrumEntropy(data)
% Input：
% data：The segment of signal
% Output：
% psdE：Power Spectrum Entropy
[pxx] = periodogram(data); 
psdE = kInformationEntopy( pxx,length(pxx) );
end

