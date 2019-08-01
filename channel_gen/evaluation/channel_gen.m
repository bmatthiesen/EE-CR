"aP:wnN1 = 2;
N2 = 1;
PL = 1;

drops = 100000;

h3 = PL*abs(1/sqrt(2)*(randn(N1, N2, drops) + 1i * randn(N1, N2, drops)));

figure
hold on;
%cdfplot(reshape(h1, [], 1));
%cdfplot(reshape(h2, [], 1));
cdfplot(reshape(h3, [], 1));
plot(0:.01:10, cdf(makedist('Rayleigh', 'b', PL/sqrt(2)), 0:.01:10));

legend('h3', 'Theoretical')
