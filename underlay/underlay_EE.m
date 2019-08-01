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

run /projects/p_mimo/cvx/cvx_startup.m 
run underlay_config

WP = str2double(getenv('SLURM_ARRAY_TASK_ID'));
savedir = getenv('UNDERLAY_HPC_SAVEDIR');
WP_SUFFIX = getenv('UNDERLAY_WP_SUFFIX');

filename = ['underlay_WPs' WP_SUFFIX '.mat'];
load(filename);

fprintf('underlay_EE: WP file = %s\tWP = %d\tSaving to %s\n', filename, WP, savedir);

h_11s = h_11s{WP};
H_12s = H_12s{WP};
h_21s = h_21s{WP};
H_22s = H_22s{WP};

drawnow('update');

for drop = 1:size(h_11s,3) %average over the channel realizations
    h_11 = h_11s(:,:,drop);
    H_12 = H_12s(:,:,drop);
    h_21 = h_21s(:,:,drop);
    H_22 = H_22s(:,:,drop);
    
    RATE = zeros(length(loading_factor),length(SNRindB));%secondary rate with rate optimization
    POWER = zeros(length(loading_factor),length(SNRindB));%secondary power with rate optimization
    EE = zeros(length(loading_factor),length(SNRindB));%secondary EE with rate optimization

    RATE_EE = zeros(length(loading_factor),length(SNRindB));%secondary rate with EE optimization
    POWER_EE = zeros(length(loading_factor),length(SNRindB));%secondary power with EE optimization
    EE_EE = zeros(length(loading_factor),length(SNRindB));%secondary EE with EE optimization
    
    for r = 1:length(loading_factor)% first loop for load factor of the primary link
        for t = 1:length(SNRindB)% second loop for SNR of the secondary link
            fprintf('[%f] drop %d  -  load factor %g/%g: %g %%\n', now, drop, r, length(loading_factor), t/length(SNRindB));
            drawnow('update');
            
            R1 = loading_factor(r) * log2(1 + P_1*norm(h_11)^2);%primary rate requirement
            I_MAX=max([0,P_1*norm(h_11)^2/(2^R1-1)-1]);%equivalent primary interference constraint
            %--------------------------------------------------------------
            %underlay spectrum sharing with rate optimization
            %CVX modelling of the case with single-user decoding at the
            %secondary receiver (Case 1)
            cvx_begin quiet;
            variable S(N_T,N_T) symmetric semidefinite;
            maximize log_det(eye(N_R)+H_22*S*H_22.'+P_1*H_12*(h_11*h_11.')*H_12.'/norm(h_11)^2)/log(2);
            h_21.'*S*h_21<=I_MAX;
            trace(S)<=P_2(t);
            cvx_end;
            
            if R1>=log2(1+P_1*norm(H_12*h_11/norm(h_11))^2)% if condition of the case is satisfied, compute the secondary rate/power/ee
                rate=cvx_optval-log2(det(eye(N_R)+P_1*H_12*(h_11*h_11.')*H_12.'/norm(h_11)^2));
                power=trace(S)+P_c;
                ee=rate/power;
            else
                %otherwise, CVX modelling of the case with successive decoding at the secondary receiver
                cvx_begin quiet;
                variable Q(N_T,N_T) symmetric semidefinite;
                maximize log_det(eye(N_R)+H_22*Q*H_22.')/log(2);
                h_21.'*Q*h_21<=I_MAX;
                trace(Q)<=P_2(t);
                cvx_end;
                if R1<=log2(det(eye(N_R)+H_22*Q*H_22.'+P_1*H_12*(h_11*h_11.')*H_12.'/norm(h_11)^2))-log2(det(eye(N_R)+H_22*Q*H_22.'))% if condition of the case is satisfied, compute the secondary rate/power/ee
                    rate=cvx_optval;
                    power=trace(Q)+P_c;
                    ee=rate/power;
                else
                    %otherwise, CVX modelling of the case with rate splitting at the secondary transmitter and successive decoding at the secondary receiver
                    cvx_begin quiet;
                    variable W1(N_T,N_T) symmetric semidefinite;
                    variable W2(N_T,N_T) symmetric semidefinite;
                    maximize log_det(eye(N_R)+H_22*(W1+W2)*H_22.'+P_1*H_12*(h_11*h_11.')*H_12.'/norm(h_11)^2)/log(2);
                    h_21.'*(W1+W2)*h_21<=I_MAX;
                    matrix_frac(sqrt(P_1)/norm(h_11)*H_12*h_11,eye(N_R)+H_22*W2*H_22.')<=2^R1-1;%matrix_frac function only accepts real channels
                    trace(W1+W2)<=P_2(t);
                    cvx_end;
                    %compute the secondary rate/power/ee
                    rate=cvx_optval-R1;
                    power=trace(W1+W2)+P_c;
                    ee=rate/power;
                end
            end
            RATE(r,t) = rate;
            POWER(r,t) = alfa .* power;
            EE(r,t) = ee ./ alfa;
            
            %--------------------------------------------------------------
            %underlay spectrum sharing with EE optimization
            if R1>=log2(1+P_1*norm(H_12*h_11/norm(h_11))^2)%the case with single-user decoding at the secondary receiver
                % Dinkelbach's method for fractional programming
                lambda=0;%initialization of parameter
                F=1;%initialization of objective value
                iter=0;%number of iterations
                while F>10^(-3) % Dinkelbach Start
                    iter=iter+1;
                    %CVX modelling
                    cvx_begin quiet
                    variable Q(N_T,N_T) symmetric semidefinite;
                    maximize log_det(eye(N_R)+H_22*Q*H_22.'+P_1*H_12*(h_11*h_11.')*H_12.'/norm(h_11)^2)/log(2)-lambda*(trace(Q)+P_c);
                    h_21.'*Q*h_21<=I_MAX;
                    %log_det(eye(N_R)+H_22*Q*H_22.'+P_1*H_12*(h_11*h_11.')*H_12.'/norm(h_11)^2)/log(2)-log2(det(eye(N_R)+P_1*H_12*(h_11*h_11.')*H_12.'/norm(h_11)^2))>=R2;
                    trace(Q)<=P_2(t);
                    cvx_end
                    F=cvx_optval-log2(det(eye(N_R)+P_1*H_12*(h_11*h_11.')*H_12.'/norm(h_11)^2));%compute updated objective value
                    lambda=(log2(det(eye(N_R)+H_22*Q*H_22.'+P_1*H_12*(h_11*h_11.')*H_12.'/norm(h_11)^2))-log2(det(eye(N_R)+P_1*H_12*(h_11*h_11.')*H_12.'/norm(h_11)^2)))/(trace(Q)+P_c);%compute updated parameter
                end    % Dinkelbach End
                %compute the secondary rate/power/ee
                rate_ee=log2(det(eye(N_R)+H_22*Q*H_22.'+P_1*H_12*(h_11*h_11.')*H_12.'/norm(h_11)^2))-log2(det(eye(N_R)+P_1*H_12*(h_11*h_11.')*H_12.'/norm(h_11)^2));
                power_ee=trace(Q)+P_c;
                ee_ee=rate_ee/power_ee;
            else
                lambda=0;
                F=1;
                iter=0;
                while F>10^(-3)
                    iter=iter+1;
                    %CVX modelling of the case with successive decoding at the secondary receiver
                    cvx_begin quiet
                    variable T(N_T,N_T) symmetric semidefinite;
                    maximize log_det(eye(N_R)+H_22*T*H_22.')/log(2)-lambda*(trace(T)+P_c);
                    h_21.'*T*h_21<=I_MAX;
                    %log_det(eye(N_R)+H_22*T*H_22.')/log(2)>=R2;
                    trace(T)<=P_2(t);
                    cvx_end
                    F=cvx_optval;
                    lambda=log2(det(eye(N_R)+H_22*T*H_22.'))/(trace(T)+P_c);
                end
                if R1<=log2(det(eye(N_R)+H_22*T*H_22.'+P_1*H_12*(h_11*h_11.')*H_12.'/norm(h_11)^2))-log2(det(eye(N_R)+H_22*T*H_22.'))% if condition of the case is satisfied, compute the secondary rate/power/ee
                    rate_ee=log2(det(eye(N_R)+H_22*T*H_22.'));
                    power_ee=trace(T)+P_c;
                    ee_ee=rate_ee/power_ee;
                else
                    lambda=0;
                    F=1;
                    iter=0;
                    while F>10^(-3)
                        iter=iter+1;
                        %CVX modelling of the case with rate splitting at the secondary transmitter and successive decoding at the secondary receiver
                        cvx_begin quiet
                        variable W1(N_T,N_T) symmetric semidefinite;
                        variable W2(N_T,N_T) symmetric semidefinite;
                        maximize log_det(eye(N_R)+H_22*(W1+W2)*H_22.'+P_1*H_12*(h_11*h_11.')*H_12.'/norm(h_11)^2)/log(2)-lambda*(trace(W1+W2)+P_c);
                        h_21.'*(W1+W2)*h_21<=I_MAX;
                        %log_det(eye(N_R)+H_22*(W1+W2)*H_22.'+P_1*H_12*(h_11*h_11.')*H_12.'/norm(h_11)^2)/log(2)-R1>=R2;
                        matrix_frac(sqrt(P_1)/norm(h_11)*H_12*h_11,eye(N_R)+H_22*W2*H_22.')<=2^R1-1;
                        trace(W1+W2)<=P_2(t);
                        cvx_end
                        F=cvx_optval-R1;
                        lambda=(log2(det(eye(N_R)+H_22*(W1+W2)*H_22.'+P_1*H_12*(h_11*h_11.')*H_12.'/norm(h_11)^2))-R1)/(trace(W1+W2)+P_c);
                    end
                    %compute the secondary rate/power/ee
                    rate_ee=log2(det(eye(N_R)+H_22*(W1+W2)*H_22.'+P_1*H_12*(h_11*h_11.')*H_12.'/norm(h_11)^2))-R1;
                    power_ee=trace(W1+W2)+P_c;
                    ee_ee=rate_ee/power_ee;
                end
            end
            RATE_EE(r,t) = rate_ee;
            POWER_EE(r,t) = power_ee .* alfa;
            EE_EE(r,t) = ee_ee ./ alfa;
        end
    end
    
    P_c = P_cin; % since power is scaled by alfa, store original P_c
    save([savedir '/' sprintf('underlay%s__WP=%d__drop=%d.mat', WP_SUFFIX, WP, drop)]);
end

exit
