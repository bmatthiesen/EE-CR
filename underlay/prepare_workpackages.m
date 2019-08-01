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

% prepares work packages for HPC jobs

clear
close all
clc;

%addpath(fileparts(pwd));

load ../channels_underlay.mat

run underlay_config

%partition_size = floor(drops / 100); % work package size roughly 100
% partition_size = floor(drops/(16*3)); % 3 jobs per cpu ???
partition_size = floor(drops / 37); % work package size roughly 37

partition = repmat(floor(drops/partition_size), 1, partition_size);
partition(1:rem(drops, partition_size)) = partition(1:rem(drops, partition_size)) + 1;

% out-of-system interference
addpath('../channel_gen');

% out-of-system base station
P_IF1_dBW = [-10 0]; % dBW
pos_IF1 = [2*p.cell_radius 0]';

% cell edge user
P_IF2_dBW = [-30 -20]; % dBW
pos_IF2 = [p.cell_radius+100 0]';

ple_IF = p.ple;


for ii = 1:length(P_IF1_dBW)
    for jj = 1:length(P_IF2_dBW)
        P_IF1 = 10^(P_IF1_dBW(ii)/10);
        P_IF2 = 10^(P_IF2_dBW(jj)/10);
        
        %% channels
        norm3 = @(E) arrayfun(@(idx) norm(E(:,:,idx)), 1:size(E,3)); % page-wise norm
        rshp = @(mat, dim) repmat(shiftdim(mat,-1),size(dim,1),size(dim,2),1);
        PL = @(d) (pathloss(d, ple_IF) .* shadowing(p.shadowing_sigma_dB, size(d))).^2;
        mkchan = @(H, N) sqrt(1./rshp(N,H)) .* H(:,:,1:drops);
        
        dist_Rx1_IF1 = norm3(p.positions.Rx1(:,:,1:drops)-repmat(pos_IF1, 1, 1, drops));
        dist_Rx1_IF2 = norm3(p.positions.Rx1(:,:,1:drops)-repmat(pos_IF2, 1, 1, drops));
        Nvar1 = Pn + P_IF1 * PL(dist_Rx1_IF1) + P_IF2 * PL(dist_Rx1_IF2);
        
        dist_Rx2_IF1 = norm3(p.positions.Rx2(:,:,1:drops)-repmat(pos_IF1, 1, 1, drops));
        dist_Rx2_IF2 = norm3(p.positions.Rx2(:,:,1:drops)-repmat(pos_IF2, 1, 1, drops));
        Nvar2 = Pn + P_IF1 * PL(dist_Rx2_IF1) + P_IF2 * PL(dist_Rx2_IF2);
        
        h_11s = mkchan(p.channels.h11, Nvar1);
        h_21s = mkchan(p.channels.h21, Nvar1);
        H_12s = mkchan(p.channels.H12, Nvar2);
        H_22s = mkchan(p.channels.H22, Nvar2);
        
        %% create work packages
        h_11s = squeeze(mat2cell(h_11s, size(h_11s,1), size(h_11s,2), partition));
        H_12s = squeeze(mat2cell(H_12s, size(H_12s,1), size(H_12s,2), partition));
        h_21s = squeeze(mat2cell(h_21s, size(h_21s,1), size(h_21s,2), partition));
        H_22s = squeeze(mat2cell(H_22s, size(H_22s,1), size(H_22s,2), partition));
        
        P_IF1 = P_IF1_dBW(ii);
        save(sprintf('underlay_WPs_PIF1=%ddBW_d1=%dm_PIF2=%ddBW_d2=%dm.mat', P_IF1_dBW(ii), pos_IF1(1), P_IF2_dBW(jj), pos_IF2(1)), 'h_11s', 'H_12s', 'h_21s', 'H_22s', 'ple_IF', 'pos_IF1', 'pos_IF2', 'P_IF1', 'P_IF2', 'dist_Rx1_IF1', 'dist_Rx1_IF2', 'Nvar1', 'dist_Rx2_IF1', 'dist_Rx2_IF2', 'Nvar2');
    end
end

fprintf('Job contains %d work packages\n', length(h_11s));
