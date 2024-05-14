function dx = Copy_of_robot_sim(x,u,p,t,data)
% Two-link robot arm Problem - Dynamics - simulation
%
% Syntax:  
%          [dx] = Dynamics(x,u,p,t,vdat)
% 
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%    vdat - structured variable containing the values of additional data used inside
%          the function%      
% Output:
%    dx - time derivative of x
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

%% States and inputs

x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);
x4 = x(:,4);
x5 = x(:,5);
x6 = x(:,6);
x7 = x(:,7);

u1 = u(:,1);
u2 = u(:,2);
%% Constants

m = 5000; %kg
Iy = 3e5;
g = 9.81;
S = 200;

k = 1.8;
cbar = 20;

Z_T = 8.*k.*(u1.^2)./m;
V_in = 0.5.*(sqrt(2*Z_T./(8*1.225*5)+x4.^2) + x4.^2);
V = sqrt(x3.^2 + x4.^2);
qbar = 0.5.*1.225.*V.^2;

%This is a simplified version, but a lot more correct than constant
C_L = -2*.1*x3./V - 1.2*x7;
C_D = -2*.012*x3./V + 0.1*x7;
C_theta = 2*.2*x3./V;

% This is even more simplified (constant)
% C_L = qbar.*S.*1.2;
% C_D = qbar.*S.*0.1;
% C_theta = qbar.*S.*cbar.*.2;

% C_L = 0;
% C_D = 0;
% C_theta = 0;


%% Forces
X_A = qbar.*S.*C_D;
% X_T = 2.*rho.*A.*V_tot_9.*V_out_9; This is correct but will try later
X_T = k.*(u2.^2)./m;
X = X_A + X_T;

Z_A = qbar.*S*C_L;
% Z_T = 2.*rho.*A.*8.*(V_tot_n.*V_out_n); This is correct will try later
Z = Z_A + Z_T;

M_A = qbar.*S.*C_theta.*cbar;
M_T = 0;
M = M_A + M_T;

%% Dynamics
dx(:,1) = x3;
dx(:,2) = x4;
dx(:,3) = (-x6.*x4 + X)./m;
dx(:,4) = (x6.*x3 + Z)./m;
dx(:,5) = x6;
dx(:,6) = M./Iy;
dx(:,7) = x4./V;