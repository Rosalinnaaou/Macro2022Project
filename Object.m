% The objective function for the PEAs described in the article 
% "Parameterized Expectations Algorithm: How to Solve for Labor Easily". 
%
% By Lilia Maliar and Serguei Maliar
% Universidad de Alicante
% February 3, 2004
% 
%-------------------------------------------------------------------

function y = object(ksiz,X)
global T k tet bt y;
y=exp(X*ksiz);
