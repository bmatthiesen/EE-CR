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

addpath(fileparts(pwd));
[~, GITHEAD] = system('git rev-parse HEAD');

savedir = 'results_testrun_24Jul';

%% config
master = 0;

drops = 1000;
P2dBW = -40:2:0;
P1dBWrange = [-30];
R1fac = [1.25 2];

problem_types = 3;

Pc = 1;
alpha = 10;
objtol = 1e-2;
feastol = 0;
MaxIter = 100;

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
%WPs = allcomb(1:drops, 1:length(R1fac), 1:length(P1range), 1:problem_types);
load(fullfile(savedir, 'WPs.mat'));
fprintf('%d work packages\n', size(WPs,1));

%% initialization
WP_init = -1*ones(1, size(WPs, 2));

%cvx_solver mosek % mosek needs jvm which makes matlab segfault when MPI is used
cvx_solver SDPT3

COM = MPI_COMM_WORLD;

TAG_WP = 1;
TAG_WORK_DONE = 2;
TAG_BYE = 3;
TAG_ERROR = -1;

[~, flag] = MPI_Initialized;
if ~flag
	MPI_Init();
end

[~, world_size] = MPI_Comm_size(COM);
[~, world_rank] = MPI_Comm_rank(COM);
[~, processor_name] = MPI_Get_processor_name;

name = [processor_name '-' num2str(world_rank)];

if world_size < 2
	% no worker
	disp('world_size must be at least 2');
	MPI_Abort(COM, 0);
end

% initialize randstream (old method due to mellum)
s = RandStream.create('mt19937ar', 'Seed', sum(100*clock)*world_rank);
RandStream.setDefaultStream(s);

%% create working directory
if world_rank == master
    if exist('savedir') && ~isempty(savedir)
        disp(['Using user supplied savedir "' savedir '"']);
    else
        savedir = ['results_' strrep(datestr(now), ' ', '-')];
    end
    
	len = length(savedir);

	MPI_Bcast(len, master, COM);
	MPI_Bcast(savedir, master, COM);
else
	len = -1;
	MPI_Bcast(len, master, COM);

	savedir = repmat(' ', 1, len);
	MPI_Bcast(savedir, master, COM);
end

[~, ~] = mkdir(savedir); % suppresses warning if folder exists

logfile = fopen([savedir '/' name '.log'], 'w');

fprintf(logfile, '%s [%f]: savedir is "%s".\n', name, now, [pwd '/' savedir]);
fprintf('%s [%f]: savedir is "%s".\n', name, now, [pwd '/' savedir]);

if strcmp(get(0, 'Diary'), 'off')
    diary([savedir '/' name '.diary']);
    fprintf(logfile, '%s [%f]: Writing diary to file %s.\n', name, now, [savedir '/' name '.diary']);
    fprintf('%s [%f]: Writing diary to file %s.\n', name, now, [savedir '/' name '.diary']);
end

%% load channels
load channels_overlay.mat
sameChannelsForR1fac = true;

if world_rank == master
    fprintf('%s [%f]: I am master!\nSelecting channels...\n', name, now);
	fprintf(logfile, '%s [%f]: I am master!\nSelecting channels...\n', name, now);
    
    make_indices = true;
    try
        load([savedir '/channel_indices.mat'])
        
        if exist('ind') && ~isempty(ind)
            make_indices = false;
        end
    catch err
        make_indices = true;
        
        if ~strcmp(err.identifier, 'MATLAB:load:couldNotReadFile')
            rethrow(err);
        end
    end

    if make_indices
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
                        fprintf('%s [%f]: Not enough good channels for R1fac = %g. Aborting...\n', name, now, R1fac(ii));
                        MPI_Abort(COM, -2);
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
        save([savedir '/channel_indices.mat'], 'ind', 'fmax');
    end
    
    % Broadcast selection
    for ii=1:length(R1fac)
        for jj=1:length(P1range)
            for snr=1:length(P2dBW)
                buf = ind{jj,ii,snr};
                info = MPI_Bcast(buf, master, COM);
                if info; error(['MPI_Recv failed (err = ' num2str(info) ') on rank ' num2str(world_rank) '.']); MPI_Abort(COM, -1); end
            end
        end
    end
else
    ind = cell(length(P1range), length(R1fac), length(P2dBW));
    for ii=1:length(R1fac)
        for jj=1:length(P1range)
            for snr=1:length(P2dBW)
                buf = zeros(drops, 1);
                info = MPI_Bcast(buf, master, COM);
                if info; error(['MPI_Recv failed (err = ' num2str(info) ') on rank ' num2str(world_rank) '.']); MPI_Abort(COM, -1); end
                ind{jj,ii,snr} = buf;
            end
        end
    end
end

h11 = cell(size(ind));
h21 = cell(size(ind));
Ht  = cell(size(ind));
H22 = cell(size(ind));

for ii=1:length(R1fac)
    for jj=1:length(P1range)
        for snr=1:length(P2dBW)
            h11{jj,ii,snr} = p.channels.h11(:,:,ind{jj,ii,snr});
            h21{jj,ii,snr} = p.channels.h21(:,:,ind{jj,ii,snr});
            Ht{jj,ii,snr}  = p.channels.Ht (:,:,ind{jj,ii,snr});
            H22{jj,ii,snr} = p.channels.H22(:,:,ind{jj,ii,snr});
        end
    end
end

NT1 = p.NT1;
NT2 = p.NT2;
NR2 = p.NR2;

clear p ind buf fmax

% selected channel per R1fac f and drop d: h11{f}(:,:,d)

%% start job processing
MPI_Barrier(COM);

if world_rank == master % master
    save([savedir '/config.mat']);
    
	waiting_nodes = world_size - 1;
	wp_cur = 1;
    wp_num = size(WPs, 1);
	buf = WP_init;
	finished = false;

	while waiting_nodes > 0
		[info, stat] = MPI_Recv(buf, MPI_ANY_SOURCE, MPI_ANY_TAG, COM);
		if info; error(['MPI_Recv failed (err = ' num2str(info) ') on rank ' num2str(world_rank) '.']); MPI_Abort(COM, 1); end

		switch stat.tag
			case TAG_WP
				if buf(1) ~= -1
					fprintf(logfile, '%s [%f]: work package %s finished by %d.\n', name, now, mat2str(buf), stat.src);
                    fprintf('%s [%f]: work package %s finished by %d.\n', name, now, mat2str(buf), stat.src);
				end

				if ~finished
					buf = WPs(wp_cur, :);
					info = MPI_Send(buf, stat.src, TAG_WP, COM);
					if info; error(['MPI_Send failed (err = ' num2str(info) ') on rank ' num2str(world_rank) '.']); MPI_Abort(COM, 2); end

					fprintf(logfile, '%s [%f]: sent work package %s to %d (%d / %d == %.2f%%).\n', name, now, mat2str(buf), stat.src, wp_cur, wp_num, wp_cur/wp_num*100);
					fprintf('%s [%f]: sent work package %s to %d (%d / %d == %.2f%%).\n', name, now, mat2str(buf), stat.src, wp_cur, wp_num, wp_cur/wp_num*100);

					wp_cur = wp_cur + 1;
					if wp_cur > wp_num
						finished = true;
					end
				else
					info = MPI_Send(buf, stat.src, TAG_WORK_DONE, COM);
					if info; error(['MPI_Send failed (err = ' num2str(info) ') on rank ' num2str(world_rank) '.']); MPI_Abort(COM, 3); end

					fprintf(logfile, '%s [%f]: suspending rank %d.\n', name, now, stat.src);
					fprintf('%s [%f]: suspending rank %d.\n', name, now, stat.src);

					waiting_nodes = waiting_nodes - 1;
				end

			case TAG_BYE
				fprintf(logfile, '%s [%f]: Rank %d left without leave :(\n', name, now, stat.src);
				fprintf('%s [%f]: Rank %d left without leave :(\n', name, now, stat.src);

				waiting_nodes = waiting_nodes - 1;

			otherwise
				info = MPI_Send(buf, stat.src, TAG_ERROR, COM);
				if info; error(['MPI_Send failed (err = ' num2str(info) ') on rank ' num2str(world_rank) '.']); MPI_Abort(COM, 4); end

				fprintf(logfile, '%s [%f]: received bad tag %d by rank %d. Suspending it.\n', name, now, stat.tag, stat.src);
				fprintf('%s [%f]: received bad tag %d by rank %d. Suspending it.\n', name, now, stat.tag, stat.src);

				waiting_nodes = waiting_nodes - 1;
		end
	end
else % slave
	wp = WP_init;

	while true
		info = MPI_Send(wp, master, TAG_WP, COM); % request wp
		if info; error(['MPI_Send failed (err = ' num2str(info) ') on rank ' num2str(world_rank) '.']); MPI_Abort(COM, 5); end

		[info stat] = MPI_Recv(wp, master, MPI_ANY_TAG, COM);
		if info; error(['MPI_Recv failed (err = ' num2str(info) ') on rank ' num2str(world_rank) '.']); MPI_Abort(COM, 6); end

		switch stat.tag
			case TAG_WP
				fprintf(logfile, '%s [%f]: Received work package %s.\n', name, now, mat2str(wp));
				fprintf('%s [%f]: Received work package %s.\n', name, now, mat2str(wp));

                drop = wp(wpidx_drop);
                job = wp(wpidx_problem_type);
                fac = R1fac(wp(wpidx_R1fac));
                P1 = P1range(wp(wpidx_P1range));
                
                R1star = fac * log2(1 + P1/Nvar * norm(h11{wp(wpidx_P1range),wp(wpidx_R1fac)}(:,:,drop))^2);
                
                p.R1star = R1star;
                p.fac = fac;
                p.drop = drop;
                p.P1dBW = P1dBWrange(wp(wpidx_P1range));
                p.P1 = P1;
                p.Nvar = Nvar;
                p.P2dBW = P2dBW;
                p.P2 = P2;
                p.alpha = alpha;
                p.Pc = Pc;
                p.objtol = objtol;
                p.error = false;
                p.processor = name;

                p.h11 = cell(1,size(h11,3));
                p.h21 = cell(1,size(h11,3));
                p.H22 = cell(1,size(h11,3));
                p.Ht = cell(1,size(h11,3));
                for snr = 1:length(P2dBW)
                    p.h11{snr} = h11{wp(wpidx_P1range),wp(wpidx_R1fac),snr}(:,:,drop);
                    p.h21{snr} = h21{wp(wpidx_P1range),wp(wpidx_R1fac),snr}(:,:,drop);
                    p.H22{snr} = H22{wp(wpidx_P1range),wp(wpidx_R1fac),snr}(:,:,drop);
                    p.Ht{snr} = Ht{wp(wpidx_P1range),wp(wpidx_R1fac),snr}(:,:,drop);
                end
                
                if job == 1
                    p.feastol = feastol;
                end

                EE = nan * ones(1, length(P2dBW));
                Rsum = nan * ones(1, length(P2dBW));
                B = nan * ones(NT2, NT2, length(P2dBW));
                
				savefile = sprintf('___R1fac=%g___P1=%gdBW___drop=%d.mat', fac, p.P1dBW, drop);
                
				switch job
					case 1
                        p.function_name = 'overlay_rank1';
                        A = nan * ones(NT2, NT2, length(P2dBW));
                        EEiter = cell(length(P2dBW),1);
                        status = cell(length(P2dBW),1);
                        a0 = cell(length(P2dBW),1);
                        
                        for snr = 1:length(P2dBW)
                            try
                                if snr > 1
                                    [A(:,:,snr), B(:,:,snr), EE(snr), a0{snr}, Rsum(snr), EEiter{snr}, status{snr}] = overlay_rank1(p.h11{snr}, p.h21{snr}, p.Ht{snr}, p.H22{snr}, p.alpha, p.R1star, p.P1, p.P2(snr), p.Pc, p.Nvar, p.objtol, p.feastol, a0{snr-1}, true, MaxIter);
                                else
                                    [A(:,:,snr), B(:,:,snr), EE(snr), a0{snr}, Rsum(snr), EEiter{snr}, status{snr}] = overlay_rank1(p.h11{snr}, p.h21{snr}, p.Ht{snr}, p.H22{snr}, p.alpha, p.R1star, p.P1, p.P2(snr), p.Pc, p.Nvar, p.objtol, p.feastol, [], true, MaxIter);
                                end
                            catch err
                                fprintf(logfile, '%s [%f]: %s (P2 = %gdBW) failed with error message "%s".\n', p.processor, now, p.function_name, p.P2dBW(snr), err.message);
                                fprintf('%s [%f]: %s (P2 = %gdBW) failed with error message "%s".\n', p.processor, now, p.function_name, p.P2dBW(snr), err.message);
                                EE(snr) = NaN;
                                Rsum(snr) = NaN;
                                p.error = true;
                            end
                        end

                        savefile = [p.function_name savefile];
                        save([savedir '/' savefile], 'EE', 'Rsum', 'A', 'B', 'p', 'EEiter', 'status');
                        clear EE Rsum A B p snr EEiter status
						

					case 2
                        p.function_name = 'overlay_KKT';
                        X = nan * ones(NT2, NT2, length(P2dBW));
                        EEiter = cell(length(P2dBW),1);
                        status = cell(length(P2dBW),1);
                        
                        for snr = 1:length(P2dBW)
                            try
                                if snr > 1
                                    [X(:,:,snr), B(:,:,snr), EE(snr), Rsum(snr), EEiter{snr}, status{snr}] = overlay_KKT(p.h11{snr}, p.h21{snr}, p.Ht{snr}, p.H22{snr}, p.alpha, p.R1star, p.P1, p.P2(snr), p.Pc, p.Nvar, p.objtol, true, MaxIter, X(:,:,snr-1));
                                else
                                    [X(:,:,snr), B(:,:,snr), EE(snr), Rsum(snr), EEiter{snr}, status{snr}] = overlay_KKT(p.h11{snr}, p.h21{snr}, p.Ht{snr}, p.H22{snr}, p.alpha, p.R1star, p.P1, p.P2(snr), p.Pc, p.Nvar, p.objtol, true, MaxIter);
                                end
                            catch err
                                fprintf(logfile, '%s [%f]: %s (P2 = %gdBW) failed with error message "%s".\n', p.processor, now, p.function_name, p.P2dBW(snr), err.message);
                                fprintf('%s [%f]: %s (P2 = %gdBW) failed with error message "%s".\n', p.processor, now, p.function_name, p.P2dBW(snr), err.message);
                                EE(snr) = NaN;
                                Rsum(snr) = NaN;
                                p.error = true;
                            end
                        end

                        savefile = [p.function_name savefile];
                        save([savedir '/' savefile], 'EE', 'Rsum', 'X', 'B', 'p', 'EEiter', 'status');
                        clear EE Rsum X B p snr EEiter status
                        
						

					case 3
                        p.function_name = 'overlay_sumratemax';
                        X = nan * ones(NT2, NT2, length(P2dBW));
                        status = cell(length(P2dBW),1);
                        
                        for snr = 1:length(P2dBW)
                            try
                                if snr > 1
                                    [X(:,:,snr), B(:,:,snr), EE(snr), Rsum(snr), status{snr}] = overlay_sumratemax(p.h11{snr}, p.h21{snr}, p.Ht{snr}, p.H22{snr}, p.alpha, p.R1star, p.P1, p.P2(snr), p.Pc, p.Nvar, p.objtol, true, MaxIter, X(:,:,snr-1));
                                else
                                    [X(:,:,snr), B(:,:,snr), EE(snr), Rsum(snr), status{snr}] = overlay_sumratemax(p.h11{snr}, p.h21{snr}, p.Ht{snr}, p.H22{snr}, p.alpha, p.R1star, p.P1, p.P2(snr), p.Pc, p.Nvar, p.objtol, true, MaxIter);
                                end
                            catch err
                                fprintf(logfile, '%s [%f]: %s (P2 = %gdBW) failed with error message "%s".\n', p.processor, now, p.function_name, p.P2dBW(snr), err.message);
                                fprintf('%s [%f]: %s (P2 = %gdBW) failed with error message "%s".\n', p.processor, now, p.function_name, p.P2dBW(snr), err.message);
                                EE(snr) = NaN;
                                Rsum(snr) = NaN;
                                p.error = true;
                            end
                        end

                        savefile = [p.function_name savefile];
                        save([savedir '/' savefile], 'EE', 'Rsum', 'X', 'B', 'p', 'status');
                        clear EE Rsum X B p snr status
				end

				fprintf(logfile, '%s [%f]: Saved results in "%s".\n', name, now, savefile);
				fprintf('%s [%f]: Saved results in "%s".\n', name, now, savefile);


			case { TAG_WORK_DONE, TAG_ERROR }
				fprintf(logfile, '%s [%f]: master signalized work done.\n', name, now);
				fprintf('%s [%f]: master signalized work done.\n', name, now);
				break; % go to BCAST

			otherwise
				fprintf(logfile, '%s [%f]: Received unknown tag. Suspending.\n', now, name);
				fprintf('%s [%f]: Received unknown tag. Suspending.\n', now, name);
				info = MPI_Send(wp, master, TAG_BYE, COM);
				if info; error(['MPI_Send failed (err = ' num2str(info) ') on rank ' num2str(world_rank) '.']); MPI_Abort(COM, 7); end
				break;
		end
	end
end


MPI_Barrier(COM);

fprintf(logfile, '%s [%f]: Exiting.\n', name, now);
fprintf('%s [%f]: Exiting.\n', name, now);
if logfile ~= 1 || logfile ~= 2 % stdout / stderr
	fclose(logfile);
end

[~, flag] = MPI_Finalized;
if ~flag
	MPI_Finalize();
end
