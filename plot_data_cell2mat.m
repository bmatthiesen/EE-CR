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

function m = plot_data_cell2mat(c)
% pad with zeros
dim = max(cell2mat(cellfun(@size, c(:), 'UniformOutput', false)));

for ii = 1:numel(c)
    for d=1:ndims(c{ii})
        dim_c = size(c{ii});
        dim_zeros = dim_c;
        dim_zeros(d) = dim(d)-dim_c(d); % pad
        c{ii} = cat(d, c{ii}, zeros(dim_zeros));
    end
end

m = cell2mat(c(:)');
end
