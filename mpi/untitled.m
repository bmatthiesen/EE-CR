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

res_path = 'res_testrun';

load(fullfile(res_path, 'config.mat'));

wpidx_drop = 1;
wpidx_problem_type = 4;
wpidx_R1fac = 2;
wpidx_P1range = 3;

names = {'overlay_rank1', 'overlay_KKT', 'overlay_sumratemax'};

for pt = 1:3
    d = dir(fullfile(res_path, [names{1} '*.mat']));
    
    for ii = 1:length(d)
        dat = load(fullfile(res_path, d(ii).name));
        
        WP = -1*ones(1, size(WPs, 2));
    
        WP(wpidx_problem_type) = pt;
        WP(wpidx_drop) = dat.p.drop;
        WP(wpidx_R1fac) = find(R1fac==dat.p.fac);
        WP(wpidx_P1range) = find(P1range==dat.p.P1);
       
        [a,idx] = ismember(WP, WPs, 'rows');
        
        assert(a);
        
        WPs(idx,:) = [];
    end
end
