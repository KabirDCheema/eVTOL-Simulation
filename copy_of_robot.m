function [dx] = copy_of_robot(x,u,p,t,data)
% Two-link robot arm Problem - Dynamics - Internal
%
% Syntax:  
%          [dx] = Dynamics(x,u,p,t,vdat)	(Dynamics Only)
%          [dx,g_eq] = Dynamics(x,u,p,t,vdat)   (Dynamics and Eqaulity Path Constraints)
%          [dx,g_neq] = Dynamics(x,u,p,t,vdat)   (Dynamics and Inqaulity Path Constraints)
%          [dx,g_eq,g_neq] = Dynamics(x,u,p,t,vdat)   (Dynamics, Equality and Ineqaulity Path Constraints)
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
%    g_eq - constraint function for equality constraints
%    g_neq - constraint function for inequality constraints
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

% x1 = x(:,1);x2 = x(:,2);x3 = x(:,3);
% u1 = u(:,1);u2 = u(:,2);
% 
% dx(:,1) = (sin(x3).*(9/4*cos(x3).*x1.^2)+2*x2.^2 + 4/3*(u1-u2) - 3/2*cos(x3).*u2 )./ (31/36 + 9/4*sin(x3).^2);
% 
% dx(:,2) = -( sin(x3).*(9/4*cos(x3).*x2.^2)+7/2*x1.^2 - 7/3*u2 + 3/2*cos(x3).*(u1-u2) )./ (31/36 + 9/4*sin(x3).^2);
% 
% dx(:,3) = x2-x1;
% 
% dx(:,4) = x1;

%% States and inputs

x_ = x(:,1);
z = x(:,2);
u_ = x(:,3);
w = x(:,4);
theta = x(:,5);
q = x(:,6);
alpha = x(:,7);

omega_n = u(:,1);
omega_9 = u(:,2);

%% Constants

m = 5000; %kg
Iy = 3e5;
g = 9.81;
S = 200;

k = 1.8;
cbar = 20;

Z_T = 8.*k.*(omega_n.^2)./m;
V_in = 0.5.*(sqrt(2*Z_T./(8*1.225*5)+w.^2) + w.^2);
V = sqrt(u_.^2 + w.^2);
qbar = 0.5.*1.225.*V.^2;

%This is a simplified version, but a lot more correct than constant
% C_L = -2*.1*u_/V - 1.2*alpha;
% C_D = -2*.012*u_/V + 0.1*alpha;
% C_theta = 2*.2*u_/V;

% This is even more simplified (constant)
% C_L = qbar.*S.*1.2;
% C_D = qbar.*S.*0.1;
% C_theta = qbar.*S.*cbar.*.2;

C_L = 0;
C_D = 0;
C_theta = 0;


%% Forces
X_A = qbar.*S.*C_D;
% X_T = 2.*rho.*A.*V_tot_9.*V_out_9; This is correct but will try later
X_T = k.*(omega_9.^2)./m;
X = X_A + X_T;

Z_A = qbar.*S*C_L;
% Z_T = 2.*rho.*A.*8.*(V_tot_n.*V_out_n); This is correct will try later
Z = Z_A + Z_T;

M_A = qbar.*S.*C_theta.*cbar;
M_T = 0;
M = M_A + M_T;

%% Dynamics
xdot = u_;
zdot = -w;

udot = (-q.*(w + V_in) + X)./m;
wdot = (q.*u_ + Z)./m;

thetadot = theta;

qdot = M./Iy;

alphadot = (u_.*wdot - u_.*w)./(u_.^2 + w.^2);

dx = [xdot zdot udot wdot thetadot qdot alphadot];