function [c,th] = check_corr(Z, maxlag, plt)
% function [c,th] = check_corr(Z, maxlag, plt)
%
% Test the auto-correlation of the residuals
%
% [c,th] = check_corr(Z, maxlag, plt)
% where:
%    Z is the transformed integral of lambda (as returned by ks_plot)
%    maxlag is the maximum lag tested (if omitted the default value is used)
%    plt, if nonzero a plot is produced (this is the default)
%
%    c returns the auto-correlation function of the residuals
%    th is the significance threshold
%
%
% Copyright (C) Luca Citi and Riccardo Barbieri, 2010-2011.
% All Rights Reserved. See LICENSE.TXT for license details.
% {lciti,barbieri}@neurostat.mit.edu
% http://users.neurostat.mit.edu/barbieri/pphrv

if nargin < 2
    maxlag = 60;
end
if nargin < 3
    plt = true;
end

small = 0.00001;
Z = min(max(Z, small), 1-small);

N = erfinv(Z .* 2 - 1); % from uniform to gaussian

Nf = fft(N - mean(N), 2^nextpow2(length(N)));
Ns = abs(Nf) .^ 2; % power spectral density

Ncorr = real(ifft(Ns)); % auto-correlation

Ncorr = Ncorr(2:maxlag+1) ./ Ncorr(1); % normalize and keep lags we are interested in

c = Ncorr;

th = 2 / sqrt(length(Z));

if plt
    figure; hold on
    plot(c, 'k.-'); 
    plot([0 maxlag], [th th], 'm:');
    plot([0 maxlag], [-th -th], 'm:');
    axis([0 maxlag -0.6 0.6]);
end
