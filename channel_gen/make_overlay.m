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

clc; clear;

%p.drops = 300000;
p.drops = 100000;

p.NT1 = 2;
p.NT2 = 2;
p.NR2 = 2;

% dist Tx1-Rx1 in meters
p.min_d11 = 250;
p.max_d11 = 750;

% percentage range of Tx1-Rx1 
p.min_tx2 = .1;
p.max_tx2 = .9;

% dist Tx2-Rx2
p.min_d22 = 10;
p.max_d22 = 100;

p.positions = generate_positions_overlay(p);

p.ple = 3.5;
p.ple_T1 = 5;
p.shadowing_sigma_dB = 6; % standard deviation of shadowing in dB 

p.channels.seed = rng;

p.channels.h11 = zeros(p.NT1, 1, p.drops);
p.channels.h21 = zeros(p.NT2, 1, p.drops);
p.channels.Ht  = zeros(p.NT2, p.NT1, p.drops);
p.channels.H22 = zeros(p.NR2, p.NT2, p.drops);

for drop = 1:p.drops
    p.channels.h11(:,:,drop) = mk_chan(norm(p.positions.Rx1(:,:,drop) - p.positions.Tx1(:,:,drop)), p.ple_T1, p.shadowing_sigma_dB, p.NT1, 1);
    p.channels.h21(:,:,drop) = mk_chan(norm(p.positions.Rx1(:,:,drop) - p.positions.Tx2(:,:,drop)), p.ple, p.shadowing_sigma_dB, p.NT2, 1);
    p.channels.Ht(:,:,drop) = mk_chan(norm(p.positions.Tx2(:,:,drop) - p.positions.Tx1(:,:,drop)), p.ple_T1, p.shadowing_sigma_dB, p.NT2, p.NT1);
    p.channels.H22(:,:,drop) = mk_chan(norm(p.positions.Rx2(:,:,drop) - p.positions.Tx2(:,:,drop)), p.ple, p.shadowing_sigma_dB, p.NR2, p.NT2);
end

save channels_overlay.mat
