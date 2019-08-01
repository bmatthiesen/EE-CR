"aP:wnfunction pos = generate_positions_overlay(p, s)
if nargin > 1
    rng(s)
end

pos.seed = rng;

pos.Tx1 = zeros(2, 1, p.drops); % just for consistency... Tx1 === [0 0]
pos.Tx2 = zeros(2, 1, p.drops);
pos.Rx1 = zeros(2, 1, p.drops);
pos.Rx2 = zeros(2, 1, p.drops);

for drop = 1:p.drops
	%% place Rx1 uniformly in x direction with distance to Tx1 between p.min_d11 and p.max_d11
	pos.Rx1(:, :, drop) = unifrnd(p.min_d11, p.max_d11) .* [1 0].';

	%% drop Tx2 on line Tx1 - Rx1
	pos.Tx2(:, :, drop) = unifrnd(p.min_tx2, p.max_tx2) * pos.Rx1(:,:,drop);


	%% drop Rx2 above Tx2
	pos.Rx2(1, :, drop) = pos.Tx2(1, :, drop);
	pos.Rx2(2, :, drop) = unifrnd(p.min_d22, p.max_d22);
end
end
