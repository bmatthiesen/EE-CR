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

%%
N = 2;          %number of antennas at the primary transmitter
N_T = 2;     %number of antennas at the secondary transmitter
N_R = 2;     %number of antennas at the secondary receiver
P_1 = 0.1;        %power of the primary transmitter
SNRindB = -30:2:10;             %power of the secondary transmitter
P_2in = 10.^(SNRindB./10);   %power of the secondary transmitter, dB to W
P_cin = 1;  %circuit power at the secondary transmitter
alfa=10;    % inverse of amplifier efficiency
P_2 = P_2in/alfa;   % Equivalent power of the secondary transmitter
P_c = P_cin/alfa;    % Equivalen circuit power at the secondary system
loading_factor = [0.5,0.75,1];     %load factor of the primary link
drops = 1000;               % number of channel realizations DROPS

%% Noise power generation
B=180*10^3;                                     % Subcarrier bandwidth in Herz
Fn_dB=3;                                          % Receiver noise figure in dB (assumed equal for all receivers)
Fn=10^(0.1*Fn_dB);                          % Receiver noise figure (assumed equal for all receivers)
No_dBm=-174;                                  % Receiver noise power spectral density in dBm
No=(1e-3)*10^(0.1*No_dBm);           % Receiver noise power spectral density
Pn=Fn*No*B;                                     % Users' noise powers
