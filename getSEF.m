function SEF = getSEF(signal,w)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
L = 256*8;                                          % Signal Length                                     % Sampling Time Interval
Fs = 256;                                              % Sampling Frequency
Fn = Fs/2;                                              % Nyquist Frequency
Fv = linspace(0, 1, fix(L/2)+1)*Fn;                     % Frequency Vector
Iv = 1:length(Fv);                                      % Index Vector
Fourier = fft(signal)/L;                              % Fourier Transform (Along Rows Of Matrix)
Fouriers = abs(Fourier);                           % Spectrum
IntSpectrum = cumtrapz(Fv, Fouriers(Iv));               % Numeric Integration
SEF = interp1(IntSpectrum, Fv, w*IntSpectrum(end), 'linear');    % Interploate To Find ‘SEF’
end

