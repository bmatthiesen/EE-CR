"aP:wnfunction pos = generate_positions_underlay(p, s)
if nargin > 1
    rng(s)
end

pos.seed = rng;

pos.Tx1 = zeros(2, 1, p.drops); % just for consistency... Tx1 === [0 0]
pos.Tx2 = zeros(2, 1, p.drops);
pos.Rx1 = zeros(2, 1, p.drops);
pos.Rx2 = zeros(2, 1, p.drops);

for drop = 1:p.drops
    %% drop Rx1 uniformly around Tx1 (between p.min_d11 and p.cell_radius)
    [pos.Rx1(1, :, drop), pos.Rx1(2, :, drop)] = unif_drop(p.min_d11, p.cell_radius);
    
    %% drop Tx2 uniformly around Tx1
    while true
        [pos.Tx2(1, :, drop), pos.Tx2(2, :, drop)] = unif_drop(p.min_d21, p.cell_radius);
        
        % ensure Tx2 has minimum distance to Rx1
        if norm(pos.Rx1(:,:,drop)-pos.Tx2(:,:,drop)) > p.rx1_radius
            break
        end
    end
    
    %% drop Rx2 uniformly around Tx2
    while true
        [pos.Rx2(1, :, drop), pos.Rx2(2, :, drop)] = unif_drop(p.min_d22, p.max_d22);
        pos.Rx2(:,:,drop) = pos.Tx2(:,:,drop) + pos.Rx2(:,:,drop);
        
        % ensure Rx2 has minimum distance to Rx1 and is within cell
        if norm(pos.Rx1(:,:,drop)-pos.Rx2(:,:,drop)) > p.rx1_radius && norm(pos.Rx2(:,:,drop)) <= p.cell_radius
            break
        end
    end
end
end

function [x,y] = unif_drop(min_d, max_d)
	r = unifrnd(min_d, max_d);
	a = unifrnd(0, 2*pi);
	
	[x,y] = pol2cart(a, r);
end
