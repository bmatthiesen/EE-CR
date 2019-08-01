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

f=n_elements/p.drops;

d=[1-f(1,:)
    1-f(1,:)-f(2,:)
1-f(1,:)-f(2,:)-f(3,:)
1-f(1,:)-f(2,:)-f(3,:)-f(4,:)
%1-f(1,:)-f(2,:)-f(3,:)-f(4,:)-f(5,:)
    ];

fprintf('inf:\t')
for jj=1:size(d,2)
    fprintf('%g\t', f(1,jj)*100)
end
fprintf('\n');

for ii=2:size(d,1)
    fprintf('%g:\t', bins(ii));
    for jj=1:size(d,2)
        fprintf('%g\t', d(ii,jj)*100);
    end
    fprintf('\n');
end
