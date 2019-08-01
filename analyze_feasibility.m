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

P1dBWrange = [-30 -20];
R1fac = [1.25 1.5 2];
P2dBW = -60:2:40;

P2 = 10.^(P2dBW/10);
P1range = 10.^(P1dBWrange/10);

N0 = -174; % dBm/Hz (noise density)
B = 180e3; % Hz (bandwidth)
F = 3; % dB (noise figure)

Nvar = 10^(N0/10-3) * B * 10^(F/10);

alpha = 10;

load channels.mat

h11 = p.channels.h11;
h21 = p.channels.h21;
Ht = p.channels.Ht;

%% make fmax
fmax = zeros(p.drops, length(P1range), length(P2dBW));

assert(all(sort(P2dBW)==P2dBW)); % we're relying on P2dBW being ordered (for align)

parfor drop = 1:p.drops
    tmp = zeros(length(P1range), length(P2dBW));
    for P1idx = 1:length(P1range)
        P1 = P1range(P1idx);
        
        for snr = 1:length(P2dBW)
            % bounds
            R1_lb = log2(1 + P1/Nvar * norm(h11(:,:,drop))^2);
            a = norm(h11(:,:,drop))^2 * P2(snr) / alpha / (Nvar*norm(h11(:,:,drop))^2 + P1 *norm(Ht(:,:,drop)*h11(:,:,drop))^2);
            A = sqrt(a) * (h21(:,:,drop)/norm(h21(:,:,drop))) * ((Ht(:,:,drop)*h11(:,:,drop))'/norm(Ht(:,:,drop)*h11(:,:,drop)));
            R1_ub = .5 * log2(1 + P1 *norm(h11(:,:,drop))^2 / Nvar + P1 /norm(h11(:,:,drop))^2 * (abs(h21(:,:,drop)'*A*Ht(:,:,drop)*h11(:,:,drop))^2)/(Nvar + h21(:,:,drop)'*(Nvar*A*A')*h21(:,:,drop)) );
            
            tmp(P1idx, snr) = R1_ub / R1_lb;
        end
    end
    fmax(drop, :, :) = tmp;
end


%% analyze
tmp = permute(fmax,[3 2 1]);

abscnt = sum(tmp >= max(R1fac), 3)

relcnt = abscnt/size(tmp,3)

fprintf('Maximum # of drops: %d\n', min(min(abscnt)))

%% plot N(P2; R1fac)

relcnt = cell(1, length(R1fac));
leg = cell(length(P1dBWrange), length(R1fac));
tmp = permute(fmax,[3 2 1]);
for R1idx = 1:length(R1fac)
    relcnt{R1idx} = sum(tmp >= R1fac(R1idx), 3) / size(tmp,3);
    
    for ii=1:length(P1dBWrange)
        leg{ii,R1idx} = sprintf('%ddBW    %g', P1dBWrange(ii), R1fac(R1idx));
    end
end

% r = [];
% for ii=1:length(R1fac)
%     r = [r ii:length(R1fac):numel(leg)];
% end

dat = horzcat(relcnt{:});
% dat = dat(:,r);
% leg = leg';

figure
plot(P2dBW, dat)
legend(reshape(leg, numel(leg), 1), 'Location', 'EastOutside')

%% save results
keys = horzcat('P2dBW', reshape(leg, 1, numel(leg)));
keys = regexprep(keys, '\s+', ' ');
cd mpi
sav('../tex/feasibility.dat', horzcat(P2dBW', dat), keys);
cd ..
