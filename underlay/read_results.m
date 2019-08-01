% Copyright (C) 2014-2016 Bho Matthiesen, Alessio Zappone, Jing Lv
% 
% This program is used in the article:
% 
% Alessio Zappone, Bho Matthiesen, and Eduard Jorswieck, "Energy Efficiency in
% MIMO Underlay and Overlay Device-to-Device Communications and Cognitive Radio
% Systems," IEEE Transactions on Signal Processing, vol. 65, no. 4, pp.
% 1026-1041 Feb. 2017, https://doi.org/10.1109/TSP.2016.2626249
% 
% 
% License:
% This program is licensed under the GPLv2 license. If you in any way use this
% code for research that results in publications, please cite our original
% article listed above.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

function [SNR, B, loading_factor, rateMax_Rsum, rateMax_EE, rateMax_Tx2Power, eeMax_Rsum, eeMax_EE, eeMax_Tx2Power] = read_results(res_path, WP_SUFFIX)
%% read dir
d = dir(fullfile(res_path, ['underlay' WP_SUFFIX '__WP*.mat']));
assert(~isempty(d));

%% load config parameters
load(fullfile(res_path, d(1).name), 'loading_factor', 'SNRindB', 'B');

% loading_factor SNRindB

WPs = load(fullfile(res_path, ['underlay_WPs' WP_SUFFIX '.mat']), 'h_11s');
drops = squeeze(cellfun(@(h) size(h, 3), WPs.h_11s));
WPs = length(WPs.h_11s);

%% pre allocate
rateMax_Rsum2 = cell(length(loading_factor), WPs);
for r = 1:length(loading_factor)
    for w = 1:WPs
        rateMax_Rsum2{r, w} = NaN * ones(drops(w), length(SNRindB));
    end
end

rateMax_EE2 = rateMax_Rsum2;
rateMax_Pt = rateMax_Rsum2;
eeMax_Rsum2 = rateMax_Rsum2;
eeMax_EE2 = rateMax_Rsum2;
eeMax_Pt = rateMax_Rsum2;

%% read results
for ii = 1:length(d)
    dat = load(fullfile(res_path, d(ii).name), 'WP', 'drop', 'RATE', 'EE', 'RATE_EE', 'EE_EE', 'POWER', 'POWER_EE', 'P_c', 'alfa');
    
    % save data
    for r = 1:length(loading_factor)
        rateMax_Rsum2{r , dat.WP}(dat.drop, :) = dat.RATE(r, :);
        rateMax_EE2{r , dat.WP}(dat.drop, :) = dat.EE(r, :);
        rateMax_Pt{r , dat.WP}(dat.drop, :) = (dat.POWER(r, :) - dat.P_c) ./ alfa; % POWER = alpha * tr() + P_c
        
        eeMax_Rsum2{r , dat.WP}(dat.drop, :) = dat.RATE_EE(r, :);
        eeMax_EE2{r , dat.WP}(dat.drop, :) = dat.EE_EE(r, :);
        eeMax_Pt{r , dat.WP}(dat.drop, :) = (dat.POWER_EE(r, :) - dat.P_c) ./ alfa; % POWER_EE = alpha * tr() + P_c
    end
end

%% join WPs + drops
rateMax_Rsum = cell(length(loading_factor), 1);
rateMax_EE = cell(length(loading_factor), 1);
rateMax_Tx2Power = cell(length(loading_factor), 1);

eeMax_Rsum = cell(length(loading_factor), 1);
eeMax_EE = cell(length(loading_factor), 1);
eeMax_Tx2Power = cell(length(loading_factor), 1);

for r = 1:length(loading_factor)
    rateMax_Rsum{r} = cat(1, rateMax_Rsum2{r,:});
    rateMax_EE{r} = cat(1, rateMax_EE2{r,:});
    rateMax_Tx2Power{r} = cat(1, rateMax_Pt{r,:});
    
    eeMax_Rsum{r} = cat(1, eeMax_Rsum2{r,:});
    eeMax_EE{r} = cat(1, eeMax_EE2{r,:});
    eeMax_Tx2Power{r} = cat(1, eeMax_Pt{r,:});
end

SNR = SNRindB;
