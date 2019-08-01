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

function sav(fn, r, k)
fieldlen = 23; % +1.16e+00

if max(cellfun(@(x) length(x), k)) > fieldlen
    error('Maximum key length is %d', fieldlen);
end

c = @(x) ceil(fieldlen/2 + length(x)/2);
k = cellfun(@(x) sprintf('  %*s%*s', c(x), strrep(x, ' ', '_'), fieldlen-c(x), ''), k, 'UniformOutput', false);

fp = fopen(fn, 'w');
fprintf(fp, '%s\n', cell2mat(k));
fclose(fp);

save(fn, 'r', '-ascii', '-append', '-double');
end
