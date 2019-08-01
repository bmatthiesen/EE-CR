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

function [EE_RA, EE_fullrank] = read_drop(drop, path)
global R1fac P1range

d = dir(fullfile(path, ['*___drop=' num2str(drop) '.mat']));

EE_fullrank = cell(length(R1fac), length(P1range));
EE_RA = cell(length(R1fac), length(P1range));

for ii = 1:length(d)
    dat = load(fullfile(path, d(ii).name));
    
    fac_idx = find(R1fac==dat.p.fac);
    P1_idx = find(P1range==dat.p.P1);
    
    if length(fac_idx) ~= 1 || length(P1_idx) ~= 1
        error('foo');
    end
    
    
    switch dat.p.function_name
        case 'overlay_rank1'
            EE_RA{fac_idx, P1_idx}(end+1, :) = dat.EE;
            
        case 'overlay_KKT'
            EE_fullrank{fac_idx, P1_idx}(end+1, :) = dat.EE;
            
        otherwise
            continue
    end
    
    if any(isnan(dat.EE))
        fprintf('Warning: EE results contains NaNs (fac=%g, P1=%g)\n', dat.p.fac, dat.p.P1);
    end
end
end
