%%  Pendulum with Spring
%% Step 1: Conceptualize the Problem
%%
% 
% <<model.jpg>>
% 
% * Model the system dynamics
% * 1 DOF $\theta$, 1 EOM $\ddot{\theta}$
% * c.m. moves in a circular path, $\hat{e_n}, \hat{e_t}$
% * assume no losses, slender bar rotating about its endpoint, planar,
% massless linear spring
% * constants of the problem:
c.m = 0.25; % slugs
c.L = 1.5; % ft, Length of bar
c.g = 32.2; % ft/s^2
%%
%
% <<StaticsProblem.jpg>>
%
c.L1 = (sqrt(243)+12)/12; % ft, distance between pins
c.Lo = 1/4; % ft, unstretched Length of spring
c.k = c.m*c.g*(c.L/2)*(sqrt(243)/18)*15/(9*c.L1); % 3.79 lb/ft
%% Step 2-3: FBD and Coordinate Frame
%%
% 
% <<FBD.jpg>>
% 
%% Step 4: Fundamental Equations
syms m L thetadot On Ot g theta L1 L k thetaddot Lo
rA = -L*[cos(theta);sin(theta)]; % position of A
rP = [-L1;0]; % position of P
rAP = rP-rA; % vector from A to P
rAPsquared = rAP.^2;
Lspring = simplify(rAPsquared(1)+rAPsquared(2))^(1/2); % Length of Spring
ehatSpring = (rP-rA)/Lspring; % spring direction
ehatn = [cos(theta);sin(theta)]; % normal direction
ehatt = [sin(theta);-cos(theta)]; % tangential direction
Fspring = k*(Lspring-Lo)*ehatSpring; % spring force vector
% spring force in normal direction
Fspringn = simplify(Fspring(1)*ehatn(1)+Fspring(2)*ehatn(2));
% spring force in tangential direction
Fspringt = simplify(Fspring(1)*ehatt(1)+Fspring(2)*ehatt(2));
%%
% $\sum \bar{F}_n
eqn(1) = m*(L/2*thetadot^2) == On - m*g*sin(theta) + Fspringn;
%%
% $\sum \bar{F}_t$
eqn(2) = m*(L/2*thetaddot) == Ot + m*g*cos(theta) + Fspringt;
%%
% $\sum M_{z_O}$
eqn(3) = (1/3*m*L^2)*thetaddot == m*g*cos(theta)*L/2 + Fspringt*L;
%% Step 5: Knowns and Unknowns
% Knowns: m (mass), L (Length of bar), Lo (unstretched Length of spring),
% L1 (distance between pin O and P), g (gravity), k (spring constant),
% $\theta$, $\dot{\theta}$ (angular position and velocity respectively,
% state variables)
%
% Unknowns: $\ddot{\theta}$ (angular acceleration of the bar, EOM), On
% (normal reaction at pin O), Ot (tangential reaction at pin O)
%
% 3 equations, 3 unknowns, solve
%% Step 6: Constraints
% N/A
%% Step 7: Solve for the EOM
x = solve(eqn,[thetaddot,On,Ot]);
%%
%
% $$\ddot{\theta} = \frac{3g}{2L}cos(\theta)-\frac{k L_1}{m L}
%   \left(1-\frac{L_o}{\sqrt{L^2+L_1^2-2LL_1 cos(\theta)}} \right)
%   sin(\theta)$$
%
%% Step 8: Solve the EOM, Solve the Problem
syms theta(t) thetadot(t)
thetaEOM = subs(x.thetaddot,{'theta','thetadot'},{theta,thetadot});
eom = odeFunction([thetadot;thetaEOM],[theta;thetadot],m,g,L,Lo,L1,k);
[T,S] = ode45(@(t,s)eom(t,s,c.m,c.g,c.L,c.Lo,c.L1,c.k),...
    linspace(0,10,1001),[pi/18,0]);
plot(T,S(:,1),'-k',[0,10],atan(9/sqrt(243))*[1,1],'-r')
xlabel('Time, sec')
ylabel('\theta, rad')
title('Figure 1: \theta vs. Time, Red line is equilibrium angle')
%%
% *Animation*
figure(2)
rA = -c.L*[cos(S(:,1)),sin(S(:,1))]; % position of A
bar = plot([0,rA(1,1)],[0,rA(1,2)],'-sk');
hold on
spring = plot([-c.L1,rA(1,1)],[0,rA(1,2)],'-*b');
axis equal
axis manual
ylim([-1.5,0])
xlabel('x, ft')
ylabel('y, ft')
for i = 2:length(T)
    set(bar,'XData',[0,rA(i,1)],'YData',[0,rA(i,2)]);
    set(spring,'XData',[-c.L1,rA(i,1)],'YData',[0,rA(i,2)]);
    drawnow
end
%% Step 9: Does it make Sense?
% Check units, magnitude, and signs
%% Appendix
% *Attributions*
% *By-hand work*