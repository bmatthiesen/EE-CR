"aP:wnfunction Xs = shadowing(shadowing_sigma_dB, sz)
    Xs = sqrt(10.^(shadowing_sigma_dB * randn(sz) ./ 10));
