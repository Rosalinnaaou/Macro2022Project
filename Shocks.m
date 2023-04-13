% The program creating series of shocks for the PEAs described in the article 
% "Parameterized Expectations Algorithm: How to Solve for Labor Easily". 
%
% By Lilia Maliar and Serguei Maliar
% Universidad de Alicante
% February 3, 2004
% 
%-------------------------------------------------------------------

clear all;
T       = 10000;
sigma   = .01;                        % Standard deviation of errors
rho     = 0.95;                       % Persistence of technology shock
tet     = zeros(T,1)+1;                 
epsi    = randn(T,1)*sigma;
for t=2:T; 
   tet(t)=tet(t-1)^rho*exp(epsi(t));
end
save shock tet; 
