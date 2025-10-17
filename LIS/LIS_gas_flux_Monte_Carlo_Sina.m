% take first 3 stations from Feb 2017 as an example, commented out for now
% n2o = [21.65, 21.01, 17.62]'; % n2o concentration nmol/kg
% n2oeq = [10.37, 10.16, 9.84]'; % n2o equilibrium concentration nmol/kg
% kw30 = [1.77, 1.95, 2.24]'; % weighted gas transfer velocity m/d
% Fn2o30 = [20.45, 21.71, 17.89]'; % 30 day weighted n2o flux in nmol/m2/d,
% I think

% take stations 6, 7, 8, 9 as an example (3 undersaturated stations)
n2o = [7.47, 7.36, 7.54, 8.76]';
n2oeq = [7.68, 7.66, 7.91, 8.39]';
kw30 = [2.49, 2.36, 1.93, 2.89]';
Fn2o30 = [-0.56, -0.74, -0.75, 1.08]';

% largest sources of uncertainty are parameterization of k, wind speed, and
% gas concentration. For Monte Carlo, assume all of these are systematic errors.
k_rsd = 0.2; % 20% relative standard deviation of gas transfer velocity (per Wanninkhof)
wind_rsd = 0.10; % approximate uncertainty in wind speed squared; suppose 5% error in winds, then multiply to 10% because it scales with wind speed squared
n2o_rsd = 0.03; % 3% relative standard deviation, estimate of uncertainty in measured concentration and equilibrium concentration

dens = 1027; % approximate density for this seawater

% create matrices to store the output
sd_F_N2O = nan(length(n2o),1);


n=1000; % run for n iterations
mc_F_N2O = nan(length(n2o),n); % create matrix to store output of monte carlo
% this loop calculates a list of possible fluxes accounting for the
% uncertainty in k, wind speed, and concentration based on monte carlo
% distribution
for i = 1:length(n2o)
    % create matrices of k, wind speed uncertainty
    mc_kw30 = kw30(i) + k_rsd.*kw30(i).*randn(n,1); % random numbers with gaussian distribution based on uncertainty in kn2o linear scaling factor (0.251 for Sc=660 wanninkof 2014 parameterization)
    mc_wind = (1 + wind_rsd.* randn(n, 1)); %  random numbers with gaussian distribution based on uncertainty in wind speed
    mc_n2o = n2o(i) + n2o_rsd .* n2o(i) .* randn(n, 1); % random numbers with gaussian distribution based on uncertainty in n2o concentration

    % flux here is k [m/d] * n2o [nmol/kg] * dens [kg/m3] * [1 m3/1000 kg]
    % = nmol/m2/d
    mc_F_N2O(i,:) = mc_kw30 .* mc_wind .* (mc_n2o - repmat(n2oeq(i),n,1)) .* dens ./1000;

    sd_F_N2O(i) = std(mc_F_N2O(i,:)); % standard deviation of fluxes

end

% get average % RSD for the all stations
rsd_F_N2O = sd_F_N2O./Fn2o30 .*100;

% output flux mean, stdev, median, percentiles 2.5, 25, 75, 97.5 for each
% station
[mean(mc_F_N2O') 
    std(mc_F_N2O') 
    median(mc_F_N2O') 
    prctile(mc_F_N2O',2.5) 
    prctile(mc_F_N2O',5)     
    prctile(mc_F_N2O',25) 
    prctile(mc_F_N2O',75) 
    prctile(mc_F_N2O',95)     
    prctile(mc_F_N2O',97.5)]