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
eqn(2) = (thetaddot) == (-m*g*sin(theta))/(m*L) - (165*10^-3)*thetadot*(abs(thetadot));

x = solve(eqn,[T,thetaddot]);

syms theta(t) thetadot(t)
thetaEOM = subs(x.thetaddot,{'theta','thetadot'},...
               {theta,thetadot});
eom = odeFunction([thetadot;thetaEOM],[theta;thetadot],g,L);

[Time,S,TE,SE,IE] = ode45(@(t,s)eom(t,s,c.g,c.L),linspace(0,100,1001),[(15*pi/180),0],options);

figure
plot(Time,S(:,1),'-k')

%x = ln(s)/t
mean_period_time = mean(TE(2:end)-TE(1:end-1));
mean_decay_rate = mean(log(SE(:,1))./TE);
disp([Time,exp(mean_decay_rate)*cos(Time*(2*pi)/mean_period_time)])
disp(mean_decay_rate)

xlabel('Time, sec')
ylabel('\theta, rad')
hold on
plot(TE,SE)
plot(Time,exp(mean_decay_rate)*cos(Time*(2*pi)/mean_period_time),'-r')
hold off

function [value isterminal direction] = event(t,s)
    % value is a function that is zero at the event
    % isterminal is 1 if you desire to terminate integration at the event
    %               0 if you desire to continue integration
    % direction defines the slope of the function value at the event
    %      1 for positive slope, -1 for negative slope, 0 for either
    % event 1: max position    

    value = s(2);
    isterminal = false;
    direction = -1;
end