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

function [ee, varargout] = overlay_fix_nans(ee, varargin)
    nans = find(isnan(ee));

    for ii=1:length(nans)
        [m,n] = ind2sub(size(ee), nans(ii));

        if n ~= 1
            ee(m,n) = ee(m,n-1);
            
            for k = 1:length(varargin)
                varargin{k}(m,n) = varargin{k}(m,n-1);
            end
        end
    end

    varargout = varargin;
end
