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

function [EE, EE_mean, EE_std, Rsum, Rsum_mean, EEiter] = read_results(name, res_path, selection_method)
global R1fac P1range P2dBW drops ind

if nargin < 3
    selection_method = 1;
end

if nargout >= 6
    iter = true;
else
    iter = false;
end

d = dir(fullfile(res_path, [name '*.mat']));

assert(~isempty(d));

% pre-allocate
EE = cell(length(R1fac), length(P1range));
for ii = 1:length(R1fac)
    for jj = 1:length(P1range)
        EE{ii,jj} = NaN * ones(drops, length(P2dBW));
    end
end
Rsum = EE;

if iter
    EEiter = EE;
end

% ordered read
for ii = 1:length(d)
    dat = load(fullfile(res_path, d(ii).name));
    
    fac_idx = find(R1fac==dat.p.fac);
    P1_idx = find(P1range==dat.p.P1);
    
    if length(fac_idx) ~= 1 || length(P1_idx) ~= 1
        error('foo');
    end
    
    EE{fac_idx, P1_idx}(dat.p.drop, :) = dat.EE;
    Rsum{fac_idx, P1_idx}(dat.p.drop, :) = dat.Rsum;
    
    if iter
        EEiter{fac_idx, P1_idx}(dat.p.drop, :) = cellfun(@(x) length(x), dat.EEiter).';
    end
end

% throw away unprocessed drops
processed_drops = cell(length(R1fac), length(P1range));
for ii = 1:length(R1fac)
    for jj = 1:length(P1range)
        processed_drops{ii,jj} = ~all(isnan(EE{ii,jj}),2);
        assert(all(processed_drops{ii,jj} == ~all(isnan(Rsum{ii,jj}),2)));
        
        EE{ii,jj} = EE{ii,jj}(processed_drops{ii,jj},:);
        Rsum{ii,jj} = Rsum{ii,jj}(processed_drops{ii,jj},:);
        
        if iter
            EEiter{ii,jj} = EEiter{ii,jj}(processed_drops{ii,jj},:);
        end
    end
end

% smooth out numerical problems
for ii = 1:length(R1fac)
    for jj = 1:length(P1range)
        [EE{ii,jj}, Rsum{ii,jj}] = overlay_fix_nans(EE{ii,jj}, Rsum{ii,jj});
    end
end

% change channel selection (this relies on the channel selection introduced into mpi/manager.m in commits a4ffc460bd0 and a5c6b9af185)
switch selection_method
    case 1
        % do nothing
        
    case 2
        % remove rows where channel becomes infeasible
        for R1idx = 1:length(R1fac)
            for P1idx = 1:length(P1range)
                l = horzcat(ind{P1idx,R1idx,:});
                l = l(processed_drops{R1idx,P1idx}, :);
                
                keep = all(bsxfun(@eq, l, l(:,1)), 2);
                EE{R1idx,P1idx} = EE{R1idx,P1idx}(keep, :);
                Rsum{R1idx,P1idx} = Rsum{R1idx,P1idx}(keep, :);
                
                if iter
                    EEiter{R1idx,P1idx} = EEiter{R1idx,P1idx}(keep, :);
                end
            end
        end
        
    case 3
        % just remove values where channel has changed
        % effectively, this means to use a good channel as long as it's feasible (looking from the high SNR regime)
        for R1idx = 1:length(R1fac)
            for P1idx = 1:length(P1range)
                l = horzcat(ind{P1idx,R1idx,:});
                l = l(processed_drops{R1idx,P1idx}, :);
                
                change = bsxfun(@ne, l, l(:,end));
                EE{R1idx,P1idx}(change) = NaN;
                Rsum{R1idx,P1idx}(change) = NaN;
                
                if iter
                    EEiter{R1idx,P1idx}(change) = NaN;
                end
            end
        end
    otherwise
        error('Unknown selection method')
end

EE_mean = cell(size(EE)); %zeros(size(SNR));
EE_std = cell(size(EE)); %zeros(size(SNR));
Rsum_mean = cell(size(Rsum));

if iter
    EEiter_mean = EE_mean;
end

for r1 = 1:length(R1fac)
    for p1 = 1:length(P1range)
        EE_mean{r1,p1} = zeros(size(P2dBW));
        Rsum_mean{r1,p1} = zeros(size(P2dBW));
        
        if iter
            EEiter_mean{r1,p1} = zeros(size(P2dBW));
        end
        
        for ii = 1:length(P2dBW)
            EE_mean{r1,p1}(ii) = mean(EE{r1,p1}( ~isnan(EE{r1,p1}(:,ii)), ii));
            EE_std{r1,p1}(ii)  =  std(EE{r1,p1}( ~isnan(EE{r1,p1}(:,ii)), ii));
            Rsum_mean{r1,p1}(ii) = mean(Rsum{r1,p1}( ~isnan(Rsum{r1,p1}(:,ii)), ii));
            
            if iter
                EEiter_mean{r1,p1}(ii) = mean(EEiter{r1,p1}( ~isnan(EEiter{r1,p1}(:,ii)), ii));
            end
        end
    end
end

if iter
    EEiter = EEiter_mean;
end
end
