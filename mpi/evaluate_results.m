% Copyright (C) 2014-2016 Bho Matthiesen, Alessio Zappone
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

res_path = 'results';
global R1fac P1range P2dBW drops ind

load(fullfile(res_path, 'config.mat'));
load(fullfile(res_path, 'channel_indices.mat'));

%% legend template
leg_temp = cell(length(R1fac), length(P1dBWrange));
for r1 = 1:length(R1fac)
    for p1 = 1:length(P1dBWrange)
        leg_temp{r1,p1} = sprintf('%gdBW   %g', P1dBWrange(p1), R1fac(r1));
    end
end

pltreshape = @(x) reshape(x, length(R1fac)*length(P1dBWrange), 1);
leg_temp = pltreshape(leg_temp);

%% overlay_rank1
overlay_rank1_legend = leg_temp;
for ii = 1:length(leg_temp)
    overlay_rank1_legend{ii} = ['Rank1  ' overlay_rank1_legend{ii}];
end

[overlay_rank1_EE, overlay_rank1_EE_mean, overlay_rank1_EE_std, overlay_rank1_Rsum, overlay_rank1_Rsum_mean, overlay_rank1_EEiter] = read_results('overlay_rank1', res_path);

overlay_rank1_plot_EE = cell2mat(pltreshape(overlay_rank1_EE_mean)).';
overlay_rank1_plot_Rsum = cell2mat(pltreshape(overlay_rank1_Rsum_mean)).';

%% overlay_KKT
overlay_KKT_legend = leg_temp;
for ii = 1:length(leg_temp)
    overlay_KKT_legend{ii} = ['KKT    ' overlay_KKT_legend{ii}];
end

[overlay_KKT_EE, overlay_KKT_EE_mean, overlay_KKT_EE_std, overlay_KKT_Rsum, overlay_KKT_Rsum_mean, overlay_KKT_EEiter] = read_results('overlay_KKT', res_path);

overlay_KKT_plot_EE = cell2mat(pltreshape(overlay_KKT_EE_mean)).';
overlay_KKT_plot_EEstd = cell2mat(pltreshape(overlay_KKT_EE_std)).';
overlay_KKT_plot_Rsum = cell2mat(pltreshape(overlay_KKT_Rsum_mean)).';

%% overlay_sumratemax
overlay_sumratemax_legend = leg_temp;
for ii = 1:length(leg_temp)
    overlay_sumratemax_legend{ii} = ['SR    ' overlay_sumratemax_legend{ii}];
end

[overlay_sumratemax_EE, overlay_sumratemax_EE_mean, overlay_sumratemax_EE_std, overlay_sumratemax_Rsum, overlay_sumratemax_Rsum_mean] = read_results('overlay_sumratemax', res_path);

overlay_sumratemax_plot_EE = cell2mat(pltreshape(overlay_sumratemax_EE_mean)).';
overlay_sumratemax_plot_EEstd = cell2mat(pltreshape(overlay_sumratemax_EE_std)).';
overlay_sumratemax_plot_Rsum = cell2mat(pltreshape(overlay_sumratemax_Rsum_mean)).';

%% plot
figure;
plot(P2dBW, overlay_rank1_plot_EE, P2dBW, overlay_KKT_plot_EE, '-.', P2dBW, overlay_sumratemax_plot_EE, '--');
legend([overlay_rank1_legend; overlay_KKT_legend; overlay_sumratemax_legend], 'Location', 'EastOutside');
title('GEE for Overlay SS');
ylabel('GEE');
xlabel('P_{2,max} [dB]');

figure;
plot(P2dBW, overlay_rank1_plot_Rsum, P2dBW, overlay_KKT_plot_Rsum, '-.', P2dBW, overlay_sumratemax_plot_Rsum, '--');
legend([overlay_rank1_legend; overlay_KKT_legend; overlay_sumratemax_legend], 'Location', 'EastOutside');
title('Sumrate for Overlay SS');
ylabel('R_{\Sigma}');
xlabel('P_{2,max} [dB]');

%% fächerplot
figure;

cnt = 0;
for ii=1:size(overlay_KKT_EE, 1)
    for jj=1:size(overlay_KKT_EE, 2)
        cnt = cnt+1;
        subplot(size(overlay_KKT_EE, 1), size(overlay_KKT_EE, 2), cnt)
        plot(P2dBW, overlay_KKT_EE{ii,jj}, 'b');
        title(sprintf('KKT  P1=%gdBW  R1fac=%g', P1dBWrange(jj), R1fac(ii)));
    end
end

%% errorbars
% P1idx = 2;
% R1facidx = 3;
% 
% figure;
% errorbar(P2dBW, overlay_KKT_EE_mean{R1facidx,P1idx}, overlay_KKT_EE_std{R1facidx,P1idx});
% title(sprintf('KKT standard deviation  P1=%gdBW  R1fac=%g', P1dBWrange(P1idx), R1fac(R1facidx)));

%% save results
bw = B/1e6;
keys = vertcat('P2dBW', overlay_rank1_legend, overlay_KKT_legend, overlay_sumratemax_legend)';
keys = regexprep(keys, '\s+', ' ');
sav('overlay_GEE.dat', horzcat(P2dBW', bw.*overlay_rank1_plot_EE, bw.*overlay_KKT_plot_EE, bw.*overlay_sumratemax_plot_EE), keys);
sav('overlay_SR.dat', horzcat(P2dBW', overlay_rank1_plot_Rsum, overlay_KKT_plot_Rsum, overlay_sumratemax_plot_Rsum), keys);
