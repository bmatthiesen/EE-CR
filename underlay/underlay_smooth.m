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

function [ee, varargout] = underlay_smooth(ee, varargin)
    nVarargs = length(varargin);

    for ii = 1:size(ee, 2)-1
        gt = ee(:,ii) > ee(:,ii+1);
        ee(gt,ii+1) = ee(gt,ii);

        for k = 1:nVarargs
            varargin{k}(gt,ii+1) = varargin{k}(gt,ii);
        end
    end
    
    varargout = varargin;
end
