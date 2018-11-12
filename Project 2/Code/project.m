%clear everything
clc
clear all
close all

% declare constants
c.m = 0.25;
c.k = 25;
c.L = 0.5;
c.I = (1/3)*c.m*c.L^2;
c.g = 9.81;

% symbolic variables
syms L theta phi k m g thetadot phidot thetaddot phiddot I

% conversion of the equations of motion from paper work into code
% note, removing the semicolon at the end of a line outputs the current
% value into the stdout
r_ab = [L*sin(theta); -L*cos(theta); 0];
r_cd = [L*sin(phi); -L*cos(phi); 0];
r_ac = [L; 0; 0];
r_bd = r_ac + r_cd - r_ab;
M_a = L/2*sin(theta)*m*g + cross(r_ab,r_bd)*k;
M_c = L/2*sin(phi)*m*g + cross(r_ab,-r_bd)*k;

% establishing relationships
eqn(1) = I*thetaddot == M_a(3);
eqn(2) = I*phiddot == M_c(3);

% solve for accelerations
x = solve(eqn,[thetaddot; phiddot]);

% construct EOM and convert into an odeFunction
syms theta(t) phi(t) thetadot(t) phidot(t)
thetaEOM = subs(x.thetaddot, {'theta', 'thetadot'}, {theta, thetadot});
phiEOM = subs(x.phiddot, {'phi', 'phidot'}, {phi, phidot});
eom = odeFunction([theta; thetaEOM; phi; phiEOM], ...
[theta; thetadot; phi; phidot], L,k,m,g,I);

[Time,S] = ode45(@(t,s)eom(t,s,c.L,c.k,c.m,c.g,c.I),linspace(0,100,1000),[pi/12 pi/12 0 0]);
disp(S)
% values seem to be 0 often, unsure on what I am doing incorrectly
figure(1)
plot(Time,S(1),'bo')



