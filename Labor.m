% The labor function for the PEAs described in the article 
% "Parameterized Expectations Algorithm: How to Solve for Labor Easily". 
%
% By Lilia Maliar and Serguei Maliar
% Universidad de Alicante
% February 3, 2004
% 
%-------------------------------------------------------------------

function [z] = labor(n,nu,alpha,am)
z=(1-n)^(-nu)*n^alpha-am;