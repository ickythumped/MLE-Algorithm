clc;
clear all; %#ok<CLALL>
close all;

%%
data = load('rr_peaks_pp7-I.txt');
R = data; % series of times of R-events [s]
n = 11;
rss = zeros(n-1, 1);
aic = zeros(n-1,1);

%%
for nparams = 2:n
    [Thetap,Mu,Kappa,L,opt] = pplikel(R, nparams);
    Var = opt.meanRR.^3 ./ Kappa; % variance of an inverse Gaussian
    Var = 1e6 * Var; % from [s^2] to [ms^2]

    mu_without_nan = Mu;
    mu_without_nan(isnan(mu_without_nan)) = [];
    nan_length = length(R) - length(mu_without_nan);
    add_zeros = zeros(1,nan_length);
    mu_final = [add_zeros mu_without_nan];
    rr_inter = [R(1); diff(R)];
    rr_inter = rr_inter(nan_length+1:end);
    display(mean(mu_without_nan))
    display(mean(rr_inter))
    [KSdistance,Z] = ks_plot(R, L, opt.delta);
    %[corr,threshold] = check_corr(Z);
    %% AIC calculation
    diff_vec = (rr_inter.' - mu_without_nan).^2;
    rss(nparams - 1) = sum(diff_vec);
    aic(nparams - 1) = 2*opt.P + length(R)*log(rss(nparams - 1)/length(R));
    disp(rss(nparams - 1))
    disp(aic(nparams - 1))

    %% plots
    t = opt.t0 + (0:length(Mu)-1) * opt.delta;
    figure; hold on
    plot(R(2:end), 1000*diff(R), 'r*')
    %plot(t, 1000*Mu, 'b+')
    plot(R, 1000*mu_final)
    legend('RR', 'First moment of IG distribution')
    xlabel('time [s]')
    ylabel('[ms]')

end
[min_aic, index_aic] = min(aic);
disp(index_aic + 1)

delta_aic = aic - min_aic; %report

exp_aic = sum(exp(-delta_aic./2));
aic_weights = exp(-delta_aic./2)./exp_aic;

