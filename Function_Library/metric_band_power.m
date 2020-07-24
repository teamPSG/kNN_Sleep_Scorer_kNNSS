function [value, S, F] = metric_band_power(band_descr, data, Fs, movingwin, mt_freq, tapers)
%Usage:
%  [BandPowerValues, Spectrogram, FrequncyAxis] = ...
%       metric_band_power(BandDescription, Data, SamplingFrequency, ...
%                                      EpochDuration, MTFreq, Tapers)
%
%Description:
%  This function calculates power in given frequency bands using multitaper
%  power estimation.
%
%Input variables:
%  BandDescription: struct, with two fields "names" and "ranges". names is
%    a size N cell of strings with a descriptive name of the frequency
%    bands, like "delta". ranges is an Nx2 double array containing the
%    lower (inclusie) and upper (exclusive) border of the band given in
%    Hertz.
%  Data: double, 1xL record of time series
%  SamplingFrequency: double, the sampling rate given in Hz
%  EpochDuration: double, the length of calculation time window given in
%    seconds
%  MTFreq: 1x2 double, frequency range of multitaper calculation given in
%    Hz.
%  Tapers: 1x2 double, a numeric vector [TW K] where TW is the
%    time-bandwidth product and K is the number of tapers to be used (less
%    than or equal to 2TW-1).
%
%Output:
%  BandPowerValues: double, 1xN+1 values of total power in a given
%    frequency band in the calculation time window, plus the last value is
%    the total power in that time window.
%  Spectrogram: double array, spectrogram in form time x frequency as
%    calculated by the mtspecgramc function of the Chronux toolbox.
%  FrequencyAxis: double array, the frequency axis associated to
%    Spectrogram.
%
%Dependencies:
%  This function uses functions from the Chronux toolbox
%  (http://chronux.org/)
%
%Author: Tamas Kiss <kiss.t@wigner.hu>

%% Parameter calculation from input
%movingwin = [2*ep_dur ep_dur];
%movingwin = [ep_dur ep_dur];
params.fpass = mt_freq;
params.tapers = tapers;
params.Fs = Fs;

%% Band power calculation
[S, ~, F] = mtspecgramc(data, movingwin, params);
value = zeros(size(S,1), length(band_descr.ranges)+1);
for bandidx = 1:length(band_descr.ranges)
    inds= F>=band_descr.ranges(bandidx,1) & F<band_descr.ranges(bandidx,2);
    value(:, bandidx) = sum(S(:,inds),2);
end
value(:, end) = sum(S,2);

FmaxI = find(F >= Fs/2, 1, 'first');
if ~isempty(FmaxI)
    F = F(1:FmaxI);
    S = S(:, 1:FmaxI);
end
if mt_freq(2) >= Fs/2
    fprintf('metric_band_power:: Warning: frequencies are only up to %0.2f Hz instead of the required %0.2f Hz\n', F(end), mt_freq(2))
end

