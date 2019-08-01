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

function [A, B, EE, a, Rsum, EEiter, status] = overlay_rank1(h11, h21, Ht, H22, alpha, R1star, P1, P2, Pc, Nvar, objtol, feastol, a0, scale, MaxIter)
% formerly known as overlay_RA
%% setup
if nargout >= 6
    saveIter = true;
    chunksize = 25;
    EEiter = zeros(1, chunksize);
else
    saveIter = false;
end

if nargin < 13 || ~(exist('a0','var') && ~isempty(a0) && ~isnan(a0))
    a0 = norm(h11)^2 * P2 / alpha / (norm(Ht*h11)^2 * P1 + Nvar*norm(h11)^2); %(27)
    
    if scale
        a0 = a0 * Nvar;
    end
    disp('overlay_rank1: init a0');
else
    disp('overlay_rank1: using supplied a0');
end

if nargin < 14
    scale = true;
else
    assert(islogical(scale));
end

if nargin < 15
    MaxIter = inf;
end

h11_orig = h11;
h21_orig = h21;
Ht_orig = Ht;
H22_orig = H22;
Nvar_orig = Nvar;
a0_orig = a0;
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

%% Initialize
p1 = Ht*h11;
np1 = norm(p1);
np1s = np1^2;

nh11s = norm(h11)^2;
nh21s = norm(h21)^2;

AMA = P1 * np1s / nh11s + Nvar; % this is NOT A*M*A' (but just a variable name... however, it's very similar to AMA')
phi = (1 + P1/Nvar * np1s/nh11s) * norm(H22*h21)^2/nh21s;

cs = 2^(2*R1star) - 1 - P1*nh11s/Nvar;
cs_feastol = 2^(2*(R1star+feastol)) - 1 - P1*nh11s/Nvar; % TODO unused

% constants for dinkelbach
dinkelbach_tol = 1e-2 * objtol;

ii = 1;

EE = -Inf;
Rsum = -Inf;
a = NaN;
B = NaN;

%% main loop
while(true)
    fprintf('\toverlay_rank1: iteration %03d: \n', ii);
    
    a_old = a;
    B_old = B;
    
    [a, B] = dinkelbach(a0, h21, H22, AMA, NR2, NT2, nh21s, phi, alpha, P1, np1s, cs, nh11s, P2, Pc, Nvar, dinkelbach_tol);
    
    if isnan(B)
        % an error occured
        A = NaN;
        EE = NaN;
        Rsum = NaN;
        status = 'failed (CVX)';
        return
    end
    
    % objective
    EE_old = EE;
    Rsum_old = Rsum;
    
    Rsum = .5 * ( log2(det(Nvar*eye(NR2) + H22 * ((h21*h21')/nh21s * a*AMA + B) * H22')) - log2(1 + a * phi) - NR2 * log2(Nvar));
    EE = Rsum / (alpha * trace(B) + alpha * a*AMA + Pc);
    
    % save intermediate results
    if saveIter
        if length(EEiter) < ii
            EEiter = [EEiter zeros(1, chunksize)];
        end
        
        EEiter(ii) = EE;
    end
    
    if EE_old > EE % EE should be increasing
        if abs(EE - EE_old) <= objtol
            fprintf('\toverlay_rank1: iteration %03d: Objective decreased and changed by less than %g. Using last value and terminating.\n', ii, objtol);
            status = 'optimal (objtol/decreased)';
        else
            fprintf('\toverlay_rank1: iteration %03d: Objective decreased. Using last value and terminating.\n', ii, objtol);
            status = 'optimal (decreased)';
        end
                
        EE = EE_old;
        Rsum = Rsum_old;
        a = a_old;
        B = B_old;
    
        break
    else
        fprintf('\bEE = %g\n', EE);
    end
    
    if abs(EE - EE_old) <= objtol
        fprintf('\toverlay_rank1: iteration %03d: Objective changed by less than %g. Terminating.\n', ii, objtol);
        status = 'optimal (objtol)';
        break
    end
    
    if ii >= MaxIter
        fprintf('\toverlay_rank1: iteration %03d: Maximum number of iterations reached. Terminating.\n', ii);
        status = 'inaccurate (maxiter)';
        break
    end
    
    a0 = a;
    ii = ii+1;
end

if saveIter
    EEiter(ii+1:end) = [];
end

A = sqrt(a) * (h21/norm(h21)) * (p1/norm(p1))';
end

function [a, B] = dinkelbach(a0, h21, H22, AMA, NR2, NT2, nh21s, phi, alpha, P1, np1s, cs, nh11s, P2, Pc, Nvar, dinkelbach_tol)
lambda = 0; % TODO check with Alessio

F = inf;
while F > dinkelbach_tol
    % cvx
    cvx_clear
    cvx_begin quiet 
        variable a nonnegative
        variable B(NT2, NT2) semidefinite

        maximize( log_det(Nvar*eye(NR2) + H22 * (a * (h21*h21')/nh21s * AMA + B) * H22') - log(1 + a0*phi) - phi * (a - a0) / (1+ a0*phi) - NR2*log(Nvar) - lambda * (alpha * trace(B) + alpha * a * AMA + Pc) )
        subject to
        a >= cs * (Nvar + h21'*B*h21) * nh11s / (nh21s * (P1 * np1s - cs*Nvar*nh11s))
        alpha * trace(B) + alpha * a * AMA <= P2
    cvx_end
    
    if strcmpi(cvx_status, 'Inaccurate/Solved')
        %fprintf(2, '\nWarning: cvx stopped with status ''Inaccurate/Solved''. Continuing.\n');
    end
    
    if ~(strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved'))
        B = NaN;
        %keyboard
        break
    end
    
    F = cvx_optval;
    
    % update lambda
    lambda = (log_det(Nvar*eye(NR2) + H22 * (a * (h21*h21')/nh21s * AMA + B) * H22') - log(1 + a0*phi) - phi * (a - a0) / (1+ a0*phi) - NR2*log(Nvar)) / (alpha * trace(B) + alpha * a * AMA + Pc);
end
end

