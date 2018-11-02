%% Simple Pendulum
%% Step 1: Conceptualize the Problem 
% 
% <<pendulum.jpg>>
% 
% * Model the system dyanamics for a simple pendulum
% * 1 Degree of Freedom $\theta \rightarrow$ 1 equation of motion
% $\ddot{\theta}$
% * Assumptions: Particle mass, no losses, planar dynamics, massless string
% * Constants
c.g = 32.2; % ft/s^2
c.m = (5/16)/c.g; % 5 oz to slugs
c.L = 3; % ft
%% Step 2 and 3: Freebody Diagram and Coordinate Frame
% Normal-Tangential Polar coordinates were chosen since the mass moves in a
% circular path relative to a fixed point. Normal points in to the center
% of the circular path and Tangential points counter-clockwise.
%
% <<FBD.jpg>>
%
%% Step 4: $\sum \bar{F} = m\bar{a}$
syms m g L theta thetadot thetaddot T
%% 
% $\sum F_n$
eqn(1) = m*(L*thetadot^2) == T - m*g*cos(theta);
%%
% $\sum F_t$
eqn(2) = m*(L*thetaddot) == -m*g*sin(theta);
%% Step 5: Knowns and Unknowns
% Knowns:
%
% * mass (m)
% * Length (L)
% * Angular Position ($\theta$)
% * Angular Velocity ($\dot{\theta}$)
% * gravity (g)
%
% Unknowns:
%
% * Tension (T)
% * Equation of Motion ($\ddot{\theta}$)
%
% 2 equations and 2 unknowns, Solve
%% Step 6: Constraints
% N/A
%% Step 7: Solve for the Equation of Motion, ($\ddot{\theta}$)
x = solve(eqn,[T,thetaddot]);
%%
% 
% $$\ddot{\theta} = -\frac{g}{L}\sin\left(\theta\right)$$
% 
%% Step 8: Solve the Equation of Motion, Solve the Problem
syms theta(t) thetadot(t)
thetaEOM = subs(x.thetaddot,{'theta','thetadot'},...
               {theta,thetadot});
eom = odeFunction([thetadot;thetaEOM],[theta;thetadot],g,L);
[Time,S] = ode45(@(t,s)eom(t,s,c.g,c.L),linspace(0,10,1001),[pi/6,0]);
%%
% Post-Process Data
figure
plot(Time,S(:,1),'-k')
xlabel('Time, sec')
ylabel('\theta, rad')
%%
% Description of motion and answers to question
% Comparisons to real-world motion
%% Step 9: Does it Make Sense?
u = symunit;
m = m*u.slug;
g = g*u.ft/u.s^2;
L = L*u.ft;
theta = 'theta';
thetadot = 'thetadot'/u.s;
thetaddot = 'thetaddot'/u.s^2;
T = T*u.lbf;
eqn = subs(eqn);
checkUnits(eqn)
EOM = subs(x.thetaddot);
%% Appendix
