% The PEA solving the standard neoclassical stochastic growth model by
% approximating the labor function on a grid and computing hours worked 
% by interpolation as described in the article "Parameterized Expectations 
% Algorithm: How to Solve for Labor Easily". 
%
% By Lilia Maliar and Serguei Maliar
% Universidad de Alicante
% February 3, 2004
% 
% The program includes five files: 
%
% 1. Main2.m 
% 2. object.m 
% 3. labor.m
% 4. shocks.m 
% 5. csolve.m (written by Christopher Sims)
%-------------------------------------------------------------------

clear all;
CPU0 = cputime;                       % Starting time

% 1. Initialize the model parameters
% ----------------------------------
T=5000; 						              % Simulation length
alpha  = 1/3;                         % Capital share
lss    = 2/3;                         % Leisure share
pkss   = 10;                          % Capital to output ratio
pcss   = 3/4;                         % Consumption to output ratio
mu      = 1.0;                        % Risk aversion of consumption
nu      = 5.0;                        % Risk aversion of leisure

% 2.1 Solve for steady state 
% -------------------------
d=(1-pcss)/pkss; % depreciation rate           
B=(1-alpha)*pkss^((1-mu)*alpha/(1-alpha))*(pcss*(1-lss))^(-mu)*lss^nu; % utility parameter
delta=1/(1-d+alpha/pkss); % discount rate
kss=((1/delta-1+d)/alpha)^(1/(alpha-1))*(1-lss); % capital ss
yss=kss/pkss; % output ss                         
css=pcss*yss; % consumption ss                         

% 2.2 Approximate the labor function
% ---------------------------------
ass=css^(-mu)*kss^alpha*(1-alpha)/B;
for j=1:100;
   am(j)=ass*(0.5+1.5*j/100);  % grid 
   N(j) = csolve('labor',(1-lss),[],0.000001,5,nu,alpha,am(j)); % value
end;

% 3. Load the shocks
% ------------------
load shock; tet=tet(1:T,1);

% 4. Initialize time series
% -------------------------
k    = zeros(T+1,1)+kss;               % Physical capital
c    = zeros(T,1);                     % Consumption
n    = zeros(T,1);                     % Leisure
e    = zeros(T-1,1);                   % Conditional expectation

% 5. Initialize the parameters of the algorithm
% ---------------------------------------------
beta    	= [log(css^(-mu)/delta); -0.001; -0.001];  % Initial coefficients  
crate    = 0.007;                                   % Speed of moving bounds
criter  	= 1e-6;            				             % Convergence criterion
update  	= .5;            				                % Updating rate 
maxiter  = 100;                                     % Max number of iterations

% 6. The Main Loop 
% -----------------             
iteration  = 0;                               % Initially, iteration is 0
dif	     = 2e-5;					             % Initally, criterion is not satified

while (dif > criter)|(hit==1)
up_bound  = kss*(2-exp(-crate*iteration));     % Upper bound
low_bound = kss*exp(-crate*iteration);         % Lower bound
hit       = 0;               

% 6.1 Given 'beta', compute the time series
% -----------------------------------------
for t=1:T;
   c(t)=(delta*exp( beta(1) + beta(2)*log(k(t)) + beta(3)*log(tet(t))))^(-1/mu);
   at = c(t)^-mu*tet(t)*k(t)^alpha*(1-alpha)/B;
   at=at*(at>=am(1))*(at<=am(100))+am(1)*(at<am(1))+am(100)*(at>am(100));
   n(t)=interp1(am,N,at,'*linear');
   k(t+1) = ((1-d)*k(t)+tet(t)*k(t)^alpha*n(t)^(1-alpha)-c(t)); 
   if k(t+1) > up_bound  
          k(t+1) = up_bound; hit=1; 
   elseif k(t+1) < low_bound
          k(t+1) = low_bound; hit=1;
   end;
end;

% 6.2 Given simulated time series, compute the expectation part
% -------------------------------------------------------------
for t = 1:T-1;  
    e(t) = c(t+1)^(-mu)*(1-d+alpha*tet(t+1)*k(t+1)^(alpha-1)*n(t+1)^(1-alpha));
end;

% 6.3 Recompute 'beta' by using NLLS regression
% ---------------------------------------------
X=[zeros(T-1,1)+1 log(k(1:T-1)) log(tet(1:T-1))];
ksi = nlinfit(X,e,'object',beta);                 

% 6.4 Update 'beta' for the next iteration 
% ----------------------------------------
iteration
dif = norm(beta-ksi)
if dif > criter;
   beta = update*ksi + (1-update)*beta;  
else;    
   break;
end;
iteration=iteration+1;
end;

% 7.0 Computational time
% ----------------------
CPU = cputime-CPU0

% 8. Plot the time series solution y, c and k 
% -------------------------------------------
time=(1:1:T);                         
subplot(3,1,1);
plot (time,k(1:T,1)), xlabel('t'), ylabel('Capital')
title('Time series solution');
subplot(3,1,2);
plot (time,c), xlabel('t'), ylabel('Consumption')
subplot(3,1,3);
plot (time,tet.*k(1:T,1).^alpha.*n(1:T,1).^(1-alpha)), xlabel('t'), ylabel('Output')