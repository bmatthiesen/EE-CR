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

clear; clc

config=load('mpi/results_13-May-2015-19:08:54/config.mat');
load('mpi/results_13-May-2015-19:08:54/channel_indices.mat');
channels=load('mpi/results_13-May-2015-19:08:54/channels.mat');
load('mpi/results_13-May-2015-19:08:54/overlay_RA___R1fac=1.5___P1=0.1___drop=1.mat');

ci = ind{config.P1range==p.P1,config.R1fac==p.fac}(p.drop);

%snr = length(p.P2dBW);

%[A, B, EE, ~, Rsum] = overlay_rank1(p.h11{snr}, p.h21{snr}, p.Ht{snr}, p.H22{snr}, p.alpha, p.R1star, p.P1, p.P2(snr), p.Pc, p.Nvar, p.objtol, p.feastol);

% drop = 100;
% Nvar = 1e-4;
% h11 = p.channels.h11(:,:,drop);
% h21 = p.channels.h21(:,:,drop);
% H22 = p.channels.H22(:,:,drop);
% Ht = p.channels.Ht(:,:,drop);
h11 = channels.p.channels.h11(:,:,ci);
h21 = channels.p.channels.h21(:,:,ci);
H22 = channels.p.channels.H22(:,:,ci);
Ht = channels.p.channels.Ht(:,:,ci);

P2 = 10^(p.SNR(end)/10);
Nvar = 2;

objtol = p.objtol; % 1e-13;
R1star = p.fac * log2(1 + p.P1/Nvar * norm(h11)^2);

a = norm(h11)^2 * P2 / p.alpha / (Nvar * norm(h11)^2 + p.P1*norm(Ht*h11)^2);
A = sqrt(a) * (h21/norm(h21)) * ((Ht*h11)'/norm(Ht*h11));
R1_ub = .5 * log2(1 + p.P1*norm(h11)^2 / Nvar + p.P1/norm(h11)^2 * (abs(h21'*A*Ht*h11)^2)/(Nvar + h21'*(Nvar*A*A')*h21) );

assert(R1_ub >= R1star);

% [X1, B1, EE1, Rsum1, EEiter1] = overlay_KKT(h11, h21, Ht, H22, p.alpha, R1star, p.P1, P2, p.Pc, Nvar, objtol, true);
% [X2, B2, EE2, Rsum2, EEiter2] = overlay_KKT(h11, h21, Ht, H22, p.alpha, R1star, p.P1, P2, p.Pc, Nvar, objtol, false);
% 
% [X1, B1, EE1, Rsum1] = overlay_sumratemax(h11, h21, Ht, H22, alpha, R1star, 1, 1e3, 1, Nvar, objtol, true);
% [X2, B2, EE2, Rsum2] = overlay_sumratemax(h11, h21, Ht, H22, alpha, R1star, 1, 1e3, 1, Nvar, objtol, false);
% 
[A1, B1, EE1, a1, Rsum1] = overlay_rank1(h11, h21, Ht, H22, p.alpha, R1star, p.P1, P2, p.Pc, Nvar, objtol, 0, true);
[A2, B2, EE2, a2, Rsum2] = overlay_rank1(h11, h21, Ht, H22, p.alpha, R1star, p.P1, P2, p.Pc, Nvar, objtol, 0, false);

fprintf(['abs_err = ' num2str(EE1 - EE2) '\nrel_err = ' num2str((EE1-EE2)/EE1*100) '%%\n'])
assert(abs(EE1-EE2)<objtol);
