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

clear

suffixes = {'_PIF1=0dBW_d1=1000m_PIF2=-20dBW_d2=600m', '_PIF1=0dBW_d1=1000m_PIF2=-30dBW_d2=600m', '_PIF1=-10dBW_d1=1000m_PIF2=-20dBW_d2=600m', '_PIF1=-10dBW_d1=1000m_PIF2=-30dBW_d2=600m'};

R = 3;

ee = cell(length(suffixes), 1);

for untitled_ii=1:length(suffixes)
    res_path = 'underlay_OOBIF-2';
    WP_SUFFIX = suffixes{untitled_ii};
    evaluate_results;
    
    ee{untitled_ii} = eeMax_EE_mean;
end

res_path = 'results_12Jul';
WP_SUFFIX = '';
evaluate_results;

close all

% for ii=1:length(suffixes)
%     figure
%     plot(SNR',ee{ii}')
%     title(suffixes{ii})
% end

eem = shiftdim(cell2mat(shiftdim(ee,-2)), 1);

xmin = 1;
xmax = find(SNR>=0,1);

m = regexp(suffixes, '_PIF1=(-?\d+)dBW_d1=1000m_PIF2=(-?\d+)dBW_d2=600m', 'tokens');
l = cellfun(@(a) sprintf('P_{IF1} = %sdBW  P_{IF2} = %sdBW', a{1}{1}, a{1}{2}), m, 'UniformOutput', false);
l = {'NO IF' l{:}};

for R = 1:3
    figure
    plot(SNR(xmin:xmax)', eeMax_EE_mean(R,xmin:xmax)', SNR(xmin:xmax)', eem(xmin:xmax,:,R));
    title(sprintf('%d%%', loading_factor(R)*100))
    legend(l, 'Location', 'northwest')
end
