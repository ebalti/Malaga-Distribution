%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main.m
%
% Created May, 2021
% Elyes Balti
% The University of Texas at Austin
%
% If you use this code or any (modified) part of it in any publication, please cite 
% the paper: E. Balti, M. Guizani, B. Hamdaoui and B. Khalfi, 
% "Aggregate Hardware Impairments Over Mixed RF/FSO Relaying Systems With Outdated CSI," 
% in IEEE Transactions on Communications, vol. 66, no. 3, pp. 1110-1123, March 2018
%
% Contact email: ebalti@utexas.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Description
% This script produces the analytical and Monte Carlo simulations of the Cumulative Distribution Function (CDF)
% of Malaga distribution. This model is used to generate the Free Space Optical (FSO) Irradiance or also termed 
% as the atmospheric turbulences.
%% Parameters 
% T: Threshold
% alpha: positive parameter related to the effective number of large-scale cells of the scattering process
% beta: natural number and it stands for the amount of fading parameter
% Mc: number of Monte Carlo iterations
% G: real variable following a gamma distribution with E[G]=1 and parameter
% m. It represents the slow fluctuation of the LOS component.
% phiA, phiB: deterministic phases of the LOS and the coupled-to-LOS scatter terms, respectively
% b0: average power of the total scatter components is denoted by 2b0 
% rho: factor expressing the amount of scattering power coupled to the LOS component.
% OMG: average power of the LOS term
% X: large scale fluctuations
% Y: small scale fluctuations
% I: optical irradiance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc
T = -20:1:20;
T = db2pow(T);
Mc = 1e6;
alpha=4.2;
beta=5;
m=4;
G = gamrnd(m,1/m,Mc,1); 
phiA = pi/2; phiB = 0; b0 = 0.596; rho = .6; OMG = 1.32;
UL = sqrt(G)*sqrt(OMG)*exp(1i*phiA);
USC = sqrt(G)*sqrt(rho*2*b0)*exp(1i*phiB);
USp = sqrt(2*b0)/sqrt(2)*(randn(Mc,1) + 1i*randn(Mc,1));
USG = sqrt(1-rho)*USp;
X = gamrnd(alpha,1/alpha,Mc,1);
Y = abs(UL + USC + USG).^2;
I = Y.*X; 
gam = 2*b0*(1-rho);
OMGP = OMG + rho*2*b0+2*sqrt(2*b0*OMG*rho)*cos(phiA-phiB);
%% Initialize the CDF vectors
CDFi = zeros(length(T),1);
CDFm = zeros(length(T),1);

for ii=1:length(T)
%% Monte Carlo    
   for kk=1:Mc
      if I(kk) < T(ii)
          CDFm(ii) = CDFm(ii) + 1;
      end      
   end
%% Analytical    
integrand = @(x)malagaPDF(alpha,beta,gam,OMGP,x);
CDFi(ii) = integral(integrand,0,T(ii));
end

CDFm = CDFm/Mc;% Averaging over the number of Monte Carlo iterations
T = pow2db(T);

figure; hold on
plot(T,CDFi,'linewidth',2)
plot(T,CDFm,'*','linewidth',2)
xlabel('Threshold (dB)')
ylabel('Cumulative Distribution Function (CDF)')
title('Malaga Distribution')
legend('Integral-Form','Monte Carlo','location','northwest')

