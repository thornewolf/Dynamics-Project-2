% * Constants
c.g = 9.81; % ms/s^2
c.m = 0.142; % kg
c.L = .5; % ft

options = odeset('Events', @event);

% Step 4: $\sum \bar{F} = m\bar{a}$
syms m g L theta thetadot thetaddot T
% $\sum F_n$
eqn(1) = m*(L*thetadot^2) == T - m*g*cos(theta);
% $\sum F_t$
eqn(2) = (thetaddot*(m*L)) == -m*g*sin(theta) - (1.65*10^-3)*thetadot*(abs(thetadot));

x = solve(eqn,[T,thetaddot]);

syms theta(t) thetadot(t)
thetaEOM = subs(x.thetaddot,{'theta','thetadot'},...
               {theta,thetadot});
eom = odeFunction([thetadot;thetaEOM],[theta;thetadot],g,L,m);

[Time,S,TE,SE,IE] = ode45(@(t,s)eom(t,s,c.g,c.L,c.m),linspace(0,100,1001),[(15*pi/180),0],options);

figure
plot(Time,S(:,1),'-k')
xlabel('Time, sec')
ylabel('\theta, rad')
hold on
title('\theta vs Time')

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [1 1 6 2]);
fig = gcf;
print('BestFitFigure','-dpdf');

avg_freq = 2*pi / (mean(diff(TE)))

function [value isterminal direction] = event(t,s)
    % value is a function that is zero at the event
    % isterminal is 1 if you desire to terminate integration at the event
    %               0 if you desire to continue integration
    % direction defines the slope of the function value at the event
    %      1 for positive slope, -1 for negative slope, 0 for either
    % event 1: max position    

    value(1) = s(2);
    isterminal(1) = false;
    direction(1) = -1;
end