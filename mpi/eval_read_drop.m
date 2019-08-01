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
sel_drop = 300;

global R1fac P1range SNR
load(fullfile(res_path, 'config.mat'));

%% legend
leg_temp = cell(length(R1fac), length(P1range));
for r1 = 1:length(R1fac)
    for p1 = 1:length(P1range)
        leg_temp{r1,p1} = sprintf('P1=%-3g  R1fac=%-3g', P1range(p1), R1fac(r1));
    end
end

pltreshape = @(x) reshape(x, length(R1fac)*length(P1range), 1);
leg_temp = pltreshape(leg_temp);

overlay_rank1_legend = leg_temp;
overlay_KKT_legend = leg_temp;
for ii = 1:length(leg_temp)
    overlay_rank1_legend{ii} = ['Rank1  ' overlay_rank1_legend{ii}];
    overlay_KKT_legend{ii} = ['KKT    ' overlay_KKT_legend{ii}];
end

%% read drop

[overlay_rank1_EE, overlay_KKT_EE] = read_drop(sel_drop, res_path);

overlay_rank1_plot_EE = cell2mat(pltreshape(overlay_rank1_EE));
overlay_KKT_plot_EE = cell2mat(pltreshape(overlay_KKT_EE));

%% plot
figure;
plot(SNR, overlay_rank1_plot_EE, SNR, overlay_KKT_plot_EE, '-.');
legend([overlay_rank1_legend; overlay_KKT_legend], 'Location', 'EastOutside');
title(['GEE for Overlay SS (only drop ' num2str(sel_drop) ')']);
ylabel('GEE');
xlabel('P_{2,max} [dB]');
