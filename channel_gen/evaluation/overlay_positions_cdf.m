"aP:wn%% make_indices (modified) 
clear; clc
logfile = 1;

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


%%
h11 = cell(size(ind,1), size(ind,2));
h21 = cell(size(ind,1), size(ind,2));
Ht  = cell(size(ind,1), size(ind,2));
H22 = cell(size(ind,1), size(ind,2));

Tx1 = cell(size(ind,1), size(ind,2));
Tx2 = cell(size(ind,1), size(ind,2));
Rx1 = cell(size(ind,1), size(ind,2));
Rx2 = cell(size(ind,1), size(ind,2));

for ii=1:length(R1fac)
    for jj=1:length(P1range)
        [~, min_snr] = min(squeeze(cellfun(@length,ind(1,1,:))));
        
        for snr = [1:min_snr-1 min_snr+1:size(ind,3)]
            assert(all(ismember(ind{1,1,1}, ind{1,1,2})))
        end
            
            h11{jj,ii} = p.channels.h11(:,:,ind{jj,ii,min_snr});
            h21{jj,ii} = p.channels.h21(:,:,ind{jj,ii,min_snr});
            Ht{jj,ii}  = p.channels.Ht (:,:,ind{jj,ii,min_snr});
            H22{jj,ii} = p.channels.H22(:,:,ind{jj,ii,min_snr});
            
            Tx1{jj,ii} = p.positions.Tx1(:,:,ind{jj,ii,min_snr});
            Tx2{jj,ii} = p.positions.Tx2(:,:,ind{jj,ii,min_snr});
            Rx1{jj,ii} = p.positions.Rx1(:,:,ind{jj,ii,min_snr});
            Rx2{jj,ii} = p.positions.Rx2(:,:,ind{jj,ii,min_snr});
    end
end

%% positions analysis
R1facidx = 1; % doesn't matter since idx are same

cd('../mpi');
savl = @sav;
cd('../channel_gen');

legend_strings = cell(length(P1range)+1,1);
keys = cell(length(P1range)+1,1);
for P1rangeidx = 1:length(P1range)
    legend_strings{P1rangeidx,1} = sprintf('P1 = %d dBW', P1dBWrange(P1rangeidx));
    keys{P1rangeidx,1} = sprintf('%ddBW_X', P1dBWrange(P1rangeidx));
    keys{P1rangeidx,2} = sprintf('%ddBW_Y', P1dBWrange(P1rangeidx));
end
legend_strings{end,1} = 'Untransformed';
keys{end,1} = 'Untransformed_X';
keys{end,2} = 'Untransformed_Y';
% keys=keys';
% keys=keys(:);
% keys=keys';

%%% cdf_Tx1Rx1
cdf_Tx1Rx1 = cell(length(P1range)+1, 1);

figure
hold on
for P1rangeidx = 1:length(P1range)
    h = cdfplot(Rx1{P1rangeidx,R1facidx}(1,1,:));
    cdf_Tx1Rx1{P1rangeidx, 1} = extract_plot_data_from_handle(h);
end

%h = plot(p.min_d11:.01:p.max_d11, cdf(makedist('uniform', 'lower', p.min_d11, 'upper', p.max_d11), p.min_d11:.01:p.max_d11));
h = cdfplot(p.positions.Rx1(1,1,:));
cdf_Tx1Rx1{end,1} = extract_plot_data_from_handle(h);

legend(legend_strings,'Location','NW')
title('CDF of distance Tx1-Rx1 (feasible channels only)')
hold off

cdf_Tx1Rx1 = plot_data_reduce(cdf_Tx1Rx1, .1);
for ii=1:numel(cdf_Tx1Rx1)
    key = strsplit(keys{ii,1}, '_');
    key = strjoin(key(1:end-1), '_');
    savl(['cdf_Tx1Rx1_' key '.dat'], cdf_Tx1Rx1{ii}, keys(ii,:));
end


%%% cdf_Tx2Rx2
cdf_Tx2Rx2 = cell(length(P1range)+1, 1);

figure
hold on
for P1rangeidx = 1:length(P1range)
    h = cdfplot(Rx2{P1rangeidx,R1facidx}(2,1,:));
    cdf_Tx2Rx2{P1rangeidx, 1} = extract_plot_data_from_handle(h);
end

%h = plot(p.min_d22:.01:p.max_d22, cdf(makedist('uniform', 'lower', p.min_d22, 'upper', p.max_d22), p.min_d22:.01:p.max_d22),'m');
h = cdfplot(p.positions.Rx2(2,1,:));
cdf_Tx2Rx2{end,1} = extract_plot_data_from_handle(h);

legend(legend_strings,'Location','NW')
title('CDF of distance Tx2-Rx2 (feasible channels only)')
hold off

cdf_Tx2Rx2 = plot_data_reduce(cdf_Tx2Rx2, .02);
for ii=1:numel(cdf_Tx2Rx2)
    key = strsplit(keys{ii,1}, '_');
    key = strjoin(key(1:end-1), '_');
    savl(['cdf_Tx2Rx2_' key '.dat'], cdf_Tx2Rx2{ii}, keys(ii,:));
end


%%% cdf_Tx2fac
cdf_Tx2fac = cell(length(P1range)+1, 1);

figure
hold on
for P1rangeidx = 1:length(P1range)
    h = cdfplot(Tx2{P1rangeidx,R1facidx}(1,1,:) ./ Rx1{P1rangeidx,R1facidx}(1,1,:));
    cdf_Tx2fac{P1rangeidx, 1} = extract_plot_data_from_handle(h);
end

%h = plot(p.min_tx2:.01:p.max_tx2, cdf(makedist('uniform', 'lower', p.min_tx2, 'upper', p.max_tx2), p.min_tx2:.01:p.max_tx2),'m');
h = cdfplot(p.positions.Tx2(1,1,:) ./ p.positions.Rx1(1,1,:));
cdf_Tx2fac{end, 1} = extract_plot_data_from_handle(h);

legend(legend_strings,'Location','NW')
title('CDF of Tx1-Tx2 / Tx1-Rx1 (feasible channels only)')
hold off

cdf_Tx2fac = plot_data_reduce(cdf_Tx2fac, .0005);
for ii=1:numel(cdf_Tx2fac)
    key = strsplit(keys{ii,1}, '_');
    key = strjoin(key(1:end-1), '_');
    savl(['cdf_Tx2fac_' key '.dat'], cdf_Tx2fac{ii}, keys(ii,:));
end

%% Statistics
disp('Keeper rate:')
for P1rangeidx = 1:length(P1range)
    fprintf('\tP1 = %d dbW:   %g%% (%d out of %d)\n', P1dBWrange(P1rangeidx), size(h11{P1rangeidx,R1facidx},3)/p.drops, size(h11{P1rangeidx,R1facidx},3), p.drops);
end
