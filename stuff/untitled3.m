clear all
close all
fr = [0:0.001:0.5];

sf_0 = 0.012;
bw =2 ;
B = 0.5*bw*sqrt(0.5*log(2));

env_radial = exp(-.5*(log(abs(fr)/sf_0).^2)/(B).^2) ;
plot(fr, env_radial)
hold on
plot(fr, 0.5*ones(size(fr)));
grid on

sf_0*(exp(-B*sqrt(2*log(2)))) 
sf_0*(exp(B*sqrt(2*log(2))))
log2(sf_0*(exp(-B*sqrt(2*log(2))))/(sf_0*(exp(B*sqrt(2*log(2))))))

%%
clear all
syms sf_0 B Bw
(solve('Bw=log2(sf_0*(exp(-B*sqrt(2*log(2))))/(sf_0*(exp(B*sqrt(2*log(2))))))','B'))
