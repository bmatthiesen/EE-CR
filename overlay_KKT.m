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

function [X, B, EE, Rsum, EEiter, status] = overlay_KKT(h11, h21, Ht, H22, alpha, R1star, P1, P2, Pc, Nvar, objtol, scale, MaxIter, X0)
% formerly known as overlay_fullrankmatrix
% implementation of algorithm sketched on page 20f (document version March 15)
%% setup
if nargout >= 5
    saveIter = true;
    chunksize = 25;
    EEiter = zeros(1, chunksize);
else
    saveIter = false;
end

if nargin < 12
    scale = true;
else
    assert(islogical(scale));
end

if nargin < 13
    MaxIter = inf;
end

h11_orig = h11;
h21_orig = h21;
Ht_orig = Ht;
H22_orig = H22;
Nvar_orig = Nvar;
if scale
    h11 = h11 / sqrt(Nvar);
    h21 = h21 / sqrt(Nvar);
    Ht = Ht / sqrt(Nvar);
    H22 = H22 / sqrt(Nvar);
    Nvar = 1;
end

%% Get dimensions
[NT2 NT1] = size(Ht);
[NR2 tmp] = size(H22);

assert(tmp==NT2);
assert(all(size(h11)==[NT1 1]));
assert(all(size(h21)==[NT2 1]));

%% Constants
dinkelbach_tol = 1e-2 * objtol; %1e-6;
% objtol = 1e-4;

nh11s = norm(h11)^2;
nHth11s = norm(Ht*h11)^2;
cs = 2^(2*R1star) - 1 - P1*nh11s/Nvar;

%% Initialize
if nargin < 14 || ~(exist('X0','var') && ~isempty(X0) && ~any(any(isnan(X0))))
    X0 = P2/(alpha*NT2) * eye(NT2);
    disp('overlay_KKT: init X0');
else
    disp('overlay_KKT: using supplied X0');
end

lambda = 0;

%% sequential convex optimization
ii = 1;
EE = -Inf;
B = NaN;
X = X0;
while (true)
    fprintf('\toverlay_KKT: iteration %03d: \n', ii);
    
    M0 = H22' * ((Nvar*eye(NR2)+ H22*X0*H22') \ H22);

    B_old = B;
    assert(all(all(X0 == X)));
    
    %% dinkelbach
    F = inf;
    failed = false;
    while F > dinkelbach_tol
        cvx_clear
        cvx_begin quiet
            variable X(NT2, NT2) semidefinite
            variable B(NT2, NT2) semidefinite

            maximize( log_det(Nvar*eye(NR2) + H22*(X+B)*H22') - log_det(Nvar*eye(NR2) + H22 * X0 * H22') - trace(M0'*(X-X0)) - lambda * alpha*trace(X + B) - lambda*Pc )
            subject to
                ((cs+1) * P1 * nHth11s) / (Nvar*nh11s + P1*nHth11s) * h21'*X*h21 >= cs * (Nvar + h21' * (X + B) * h21)
                alpha * trace(X + B) <= P2
        cvx_end

    %     if strcmpi(cvx_status, 'Inaccurate/Solved')
    %         fprintf(2, '\nWarning: cvx stopped with status ''Inaccurate/Solved''. Continuing.\n');
    %     end

        if ~(strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved'))
            failed = true;
            break
        end

        F = cvx_optval;
        lambda = (log_det(Nvar*eye(NR2) + H22*(X+B)*H22') - log_det(Nvar*eye(NR2) + H22 * X0 * H22') - trace(M0'*(X-X0))) / (alpha*trace(X + B) + Pc);
    end
    
    if failed
        EE = NaN;
        Rsum = NaN;
        status = 'failed (CVX)';
        return
    end
    
    %% calculate objective
    EE_old = EE;
    
    EE = .5 * log2(det(eye(NR2) + ((Nvar_orig*eye(NR2) + H22_orig*X*H22_orig') \ H22_orig*B*H22_orig'))) / (alpha * trace(X + B) + Pc);
    
    %% save intermediate results
    if saveIter
        if length(EEiter) < ii
            EEiter = [EEiter zeros(1, chunksize)];
        end
        
        EEiter(ii) = EE;
    end
    
    %% status output / check that objective is increasing
    if EE_old > EE % EE should be increasing
        if abs(EE - EE_old) <= objtol
            fprintf('\toverlay_KKT: iteration %03d: Objective decreased and changed by less than %g. Using last value and terminating.\n', ii, objtol);
            status = 'optimal (objtol/decreased)';
        else
            fprintf('\toverlay_KKT: iteration %03d: Objective decreased. Using last value and terminating.\n', ii, objtol);
            status = 'optimal (decreased)';
        end
                
        EE = EE_old;
        X = X0;
        B = B_old;
    
        break
    else
        fprintf('\bEE = %g\n', EE);
    end
    
    %% check termination criterion
    if abs(EE - EE_old) <= objtol
        fprintf('\toverlay_KKT: iteration %03d: Objective changed by less than %g. Terminating.\n', ii, objtol);
        status = 'optimal (objtol)';
        break
    end
    
    if ii >= MaxIter
        fprintf('\toverlay_KKT: iteration %03d: Maximum number of iterations reached. Terminating.\n', ii);
        status = 'inaccurate (maxiter)';
        break
    end
    
    %% update variables for next loop
    X0 = X;
    ii = ii+1;
end

if saveIter
    EEiter(ii+1:end) = [];
end

Rsum = .5 * log2(det(eye(NR2) + ((Nvar_orig*eye(NR2) + H22_orig*X*H22_orig') \ H22_orig*B*H22_orig')));
end
