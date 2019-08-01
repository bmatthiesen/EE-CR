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

 res_path = 'underlay_OOBIF-2';
 WP_SUFFIX = '_PIF1=-10dBW_d1=1000m_PIF2=-30dBW_d2=600m';
res_path = 'results_12Jul';
WP_SUFFIX = '';

[SNR, B, loading_factor, rateMax_Rsum, rateMax_EE, rateMax_Tx2Power, eeMax_Rsum, eeMax_EE, eeMax_Tx2Power] = read_results(res_path, WP_SUFFIX);

% smooth
for ii = 1:length(eeMax_EE)
    [eeMax_EE{ii}, eeMax_Rsum{ii}, eeMax_Tx2Power{ii}] = underlay_smooth(eeMax_EE{ii}, eeMax_Rsum{ii}, eeMax_Tx2Power{ii});
end

% calc mean
rateMax_Rsum_mean = cell2mat(cellfun(@underlay_mean, rateMax_Rsum, 'UniformOutput', false));
rateMax_EE_mean = B/10^6 .* cell2mat(cellfun(@underlay_mean, rateMax_EE, 'UniformOutput', false));
rateMax_Tx2Power_mean = 10*log10(cell2mat(cellfun(@underlay_mean, rateMax_Tx2Power, 'UniformOutput', false)));

eeMax_Rsum_mean = cell2mat(cellfun(@underlay_mean, eeMax_Rsum, 'UniformOutput', false));
eeMax_EE_mean = B/10^6 .* cell2mat(cellfun(@underlay_mean, eeMax_EE, 'UniformOutput', false));
eeMax_Tx2Power_mean = 10*log10(cell2mat(cellfun(@underlay_mean, eeMax_Tx2Power, 'UniformOutput', false)));

% lift 100% curves
eeMax_Rsum_mean(3,:) = eeMax_Rsum_mean(3,:) + max(rateMax_Rsum_mean(3,1:8) - eeMax_Rsum_mean(3,1:8));
eeMax_EE_mean(3,:) = eeMax_EE_mean(3,:) + max(rateMax_EE_mean(3,1:8) - eeMax_EE_mean(3,1:8));

% legend
keys = cellfun(@(x) sprintf('_%g',x), num2cell(loading_factor), 'UniformOutput', false);
keys = ['SNR', cellfun(@(x) ['EE' x], keys, 'UniformOutput', false),  cellfun(@(x) ['Rate' x], keys, 'UniformOutput', false)];

% plot
pltrng = 1:length(SNR); % 1:16

figure
plot(SNR(pltrng), eeMax_Rsum_mean(:,pltrng)', SNR(pltrng), rateMax_Rsum_mean(:,pltrng)', '--')
legend(strrep(keys(2:end),'_',' '), 'Location', 'northwest')
xlabel('P_2 [dBW]')
ylabel('R_2 [bit/s/Hz]')
title('Average R_2')

figure
plot(SNR(pltrng), eeMax_EE_mean(:,pltrng)', SNR(pltrng), rateMax_EE_mean(:,pltrng)', '--')
legend(strrep(keys(2:end),'_',' '), 'Location', 'northwest')
xlabel('P_2 [dBW]')
ylabel('EE_2 [Mbit/J]')
title('Average EE_2')

figure
plot(SNR(pltrng), eeMax_Tx2Power_mean(:,pltrng)', SNR(pltrng), rateMax_Tx2Power_mean(:,pltrng)', '--')
legend(strrep(keys(2:end),'_',' '), 'Location', 'northwest')
xlabel('P_2 [dBW]')
ylabel('tr(K_{2,1} + K_{2,2}) [dBW]')
title('Average transmit power of Tx2')

%% save
cd ../mpi
savp = @sav;
cd ../underlay

savp(['underlay_EE' WP_SUFFIX '.dat'], [SNR' eeMax_EE_mean' rateMax_EE_mean'], keys);
savp(['underlay_Rate' WP_SUFFIX '.dat'], [SNR' eeMax_Rsum_mean' rateMax_Rsum_mean'], keys);
savp(['underlay_Power' WP_SUFFIX '.dat'], [SNR' eeMax_Tx2Power_mean' rateMax_Tx2Power_mean'], keys);
