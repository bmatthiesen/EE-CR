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

% this is a rewrite of the index selection for mpi/manager.m
% it was used for testing. the new code is incorporated in mpi/manager.m now

clear; clc
logfile = 1;

drops = 10000;
P2dBW = -20:2:10;
P1dBWrange = [-30 -20];
R1fac = [1.25 1.5 2];

problem_types = 3;

Pc = 1;
alpha = 10;
objtol = 1e-4;
feastol = 0;

N0 = -174; % dBm/Hz (noise density)
B = 180e3; % Hz (bandwidth)
F = 3; % dB (noise figure)

Nvar = 10^(N0/10-3) * B * 10^(F/10);

P2 = 10.^(P2dBW/10);
P1range = 10.^(P1dBWrange/10);

wpidx_drop = 1;
wpidx_problem_type = 4;
wpidx_R1fac = 2;
wpidx_P1range = 3;

name = 'foo';




load channels_overlay.mat
sameChannelsForR1fac = true;

h11 = p.channels.h11;
h21 = p.channels.h21;
Ht = p.channels.Ht;

fmax = zeros(p.drops, length(P1range), length(P2dBW));

assert(all(sort(P2dBW)==P2dBW)); % we're relying on P2dBW being ordered (for align)

for drop = 1:p.drops
    for P1idx = 1:length(P1range)
        P1 = P1range(P1idx);
        
        for snr = 1:length(P2dBW)
            % bounds
            R1_lb = log2(1 + P1/Nvar * norm(h11(:,:,drop))^2);
            a = norm(h11(:,:,drop))^2 * P2(snr) / alpha / (Nvar*norm(h11(:,:,drop))^2 + P1 *norm(Ht(:,:,drop)*h11(:,:,drop))^2);
            A = sqrt(a) * (h21(:,:,drop)/norm(h21(:,:,drop))) * ((Ht(:,:,drop)*h11(:,:,drop))'/norm(Ht(:,:,drop)*h11(:,:,drop)));
            R1_ub = .5 * log2(1 + P1 *norm(h11(:,:,drop))^2 / Nvar + P1 /norm(h11(:,:,drop))^2 * (abs(h21(:,:,drop)'*A*Ht(:,:,drop)*h11(:,:,drop))^2)/(Nvar + h21(:,:,drop)'*(Nvar*A*A')*h21(:,:,drop)) );
            
            fmax(drop, P1idx, snr) = R1_ub / R1_lb;
        end
    end
end

% select channels
if sameChannelsForR1fac
    [~, R1facrange] = max(R1fac);
else
    R1facrange = 1:length(R1fac);
end

ind = cell(length(P1range), length(R1fac), length(P2dBW));
for ii=R1facrange
    for jj=1:length(P1range)
        % select
        for snr=1:length(P2dBW)
            ind{jj,ii,snr} = find(fmax(:,jj,snr)>=R1fac(ii));
            
            if length(ind{jj,ii,snr}) < drops
                fprintf(logfile, '%s [%f]: Not enough good channels for R1fac = %g. Aborting...\n', name, now, R1fac(ii));
                %             MPI_Abort(COM, -2);
            end
            
            ind{jj,ii,snr} = ind{jj,ii,snr}(1:drops); % select enough channels
        end
        
        % align
        for snr = size(ind,3):-1:2
            new=find(~ismember(ind{jj,ii,snr-1}, ind{jj,ii,snr}),1); % index of first new index in ind{jj,ii,snr-1}
            c = ind{jj,ii,snr}; % use ind{jj,ii,snr} as base for updated ind{jj,ii,snr-1}
            c(~ismember(ind{jj,ii,snr}, ind{jj,ii,snr-1})) = ind{jj,ii,snr-1}(new:end); % replace non-feasible channels
            
            ind{jj,ii,snr-1} = c;
        end
    end
end

if sameChannelsForR1fac
    range = 1:length(R1fac);
    range(range==R1facrange) = [];
    
    for idx=range
        for ii=1:size(ind,1)
            for jj=1:size(ind,3)
                ind{ii,idx,jj} = ind{ii,R1facrange,jj};
            end
        end
    end
end

% save selection
% save([savedir '/channel_indices.mat'], 'ind');

%% validate
for jj=1:length(P1range)
    for snr=1:length(P2dBW)
        for ii=2:length(R1fac)
            res = NaN * ones(size(ind{jj,ii,snr}));
            
            for idx = 1:length(ind{jj,ii,snr})
                res(idx) = any(ind{jj,ii,snr}(idx)==ind{jj,ii-1,snr});
            end
            
            tmp = find(res==0,1);
            if isempty(tmp)
                tmp = length(res)+1;
            else
                if ~(ind{jj,ii,snr}(tmp)>ind{jj,ii-1,snr}(end))
                    fprintf('Assertion ii 2 failed: ii=%d\tjj=%d\tsnr=%d\n', ii, jj, snr);
                end
            end
            
            if ~all(find(res==1)<tmp)
                fprintf('Assertion ii failed: ii=%d\tjj=%d\tsnr=%d\n', ii, jj, snr);
            end
        end
    end
end

disp('make_indices: done')
