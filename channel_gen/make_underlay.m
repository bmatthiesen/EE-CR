"aP:wnclc; clear;

%p.drops = 300000;
p.drops = 100000;

p.NT1 = 2;
p.NT2 = 2;
p.NR2 = 2;

p.cell_radius = 500;
p.min_d11 = 50;
p.min_d21 = 0;
p.rx1_radius = 10;
p.min_d22 = 10;
p.max_d22 = 100;

p.positions = generate_positions_underlay(p);

p.ple = 3.5;
p.shadowing_sigma_dB = 6; % standard deviation of shadowing in dB 

p.channels.seed = rng;

p.channels.h11 = zeros(p.NT1, 1, p.drops);
p.channels.h21 = zeros(p.NT2, 1, p.drops);
p.channels.H12  = zeros(p.NR2, p.NT1, p.drops);
p.channels.H22 = zeros(p.NR2, p.NT2, p.drops);

for drop = 1:p.drops
    p.channels.h11(:,:,drop) = mk_chan(norm(p.positions.Rx1(:,:,drop) - p.positions.Tx1(:,:,drop)), p.ple, p.shadowing_sigma_dB, p.NT1, 1);
    p.channels.h21(:,:,drop) = mk_chan(norm(p.positions.Rx1(:,:,drop) - p.positions.Tx2(:,:,drop)), p.ple, p.shadowing_sigma_dB, p.NT2, 1);
    p.channels.H12(:,:,drop) = mk_chan(norm(p.positions.Rx2(:,:,drop) - p.positions.Tx1(:,:,drop)), p.ple, p.shadowing_sigma_dB, p.NR2, p.NT1);
    p.channels.H22(:,:,drop) = mk_chan(norm(p.positions.Rx2(:,:,drop) - p.positions.Tx2(:,:,drop)), p.ple, p.shadowing_sigma_dB, p.NR2, p.NT2);
end

save channels_underlay.mat
