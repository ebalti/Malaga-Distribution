%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% malagaPDF.m
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
% This script returns the analytical expression of the Probability Density
% Function (PDF) of the Malaga distribution
%% Parameters
% x: argument
% alpha: positive parameter related to the effective number of large-scale cells of the scattering process
% beta: natural number and it stands for the amount of fading parameter
% gam: 2b0(1 − ρ)
% omgp: (Omega Prime) Ω+ρ2b0 + 2 √2b0Ωρ cos(φA − φB).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = malagaPDF(alpha,beta,gam,omgp,x)
A = 2*alpha^(alpha/2)/(gam^(1+alpha/2)*gamma(alpha))*(gam*beta/(gam*beta+omgp))^(beta+alpha/2);
S=0;
for k=1:beta
ak = nchoosek(beta-1,k-1)*(gam*beta+omgp)^(1-k/2)/factorial(k-1)*(omgp/gam)^(k-1)*(alpha/beta)^(k/2);
S = S + ak*x.^( (alpha+k)/2-1) .*besselk(alpha-k,2*sqrt(alpha*beta*x./(gam*beta+omgp)));
end
output = A*S;
end
