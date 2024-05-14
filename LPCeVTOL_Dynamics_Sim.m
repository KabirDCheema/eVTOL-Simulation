function dx = LPCeVTOL_Dynamics_Sim(x,u,p,t,vdat)
%Van der Pal Oscillator problem Dynamics - simulation
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
%
%------------- BEGIN CODE --------------

m = 5000; %kg
Ixz = 10;
Ix = 2;
Iz = 50;
Iy = 25;
g = 9.81;
S = 200;
CLalpha = .1;
CL1 = .1;
CD1 = .012;
CYbeta = .2;
Ctheta1 = .2;
Cthetaalpha = .1;
k = 1.8;
cbar = 20;

u_ = x(:,4);
v = x(:,5);
w = x(:,6);
phi = x(:,7);
theta = x(:,8);
psi = x(:,9);
p = x(:,10);
q = x(:,11);
r = x(:,12);
alpha = x(:,13);
beta = x(:,14);

w_per = u(:,1); % omega for each prop should be the same
w9 = u(:,2);
%delR = u(:,3);
%delA = u(:,4);
%delE = u(:,5);
ZT = 8.*k.*(w_per.^2)./m;
Vin = 0.5.*(sqrt(2*ZT./(8*1.225*5)+w.^2) + w.^2);
V = sqrt(u.^2 + v.^2 + w.^2);
qbar = 0.5.*1.225.*V.^2;
% 
xdot = u_;
ydot = v;
zdot = w;
%xdot = u + (phi.*theta - psi).*v + (theta - phi.*psi).*w;
% ydot = u.*psi + (phi.*theta.*psi - 1).*v + (theta.*psi - phi.*psi).*w;
%zdot = -u.*theta + phi.*v + w;
% 
XA = qbar.*S.*(-2.*CD1.*u./V + CL1.*alpha);
XT = k.*(w9.^2)./m;
udot = -g.*sin(theta) - q.*(w + Vin) + r.*v + XT + XA./m;
% 
% YA = qbar.*S.*(CYbeta.*beta);
% vdot = g.*cos(theta).*sin(phi) - r.*u + p.*(w + Vin) + YA./m;
vdot = 0;
ZA = qbar.*S.*(2.*CL1.*u./V - CLalpha.*alpha);
wdot = g.*cos(theta).*cos(phi) - p.*v - q.*u - ZA./m - ZT; % Did not include V_in_dot
% 
% phidot = p + q.*phi.*theta + r.*theta;
phidot = 0;
thetadot = q - r.*phi;
psidot = 0;
% psidot = q.*phi + r;
% 
M = 2.*qbar.*S.*cbar.(2.*Ctheta1.*(u + Vin)./V + Cthetaalpha.*alpha);
% 
% pdot = 0;
qdot = (M + Ixz.*(r.^2 - p.^2)+ (Iz - Ix).*r.*p)./Iy;
rdot = 0;
% 
alphadot = (u_.*wdot - u_.*w)./(u_.^2 + w.^2);
% betadot = (vdot.*(u_.^2 + w.^2) - v.*(u_.^2 +w.*wdot))/(udot.^2 * sqrt(u_.^2 + w.^2));
betadot=0;



dx = [xdot ydot zdot udot vdot wdot phidot thetadot psidot pdot qdot rdot alphadot betadot];
%------------- END OF CODE ----------------