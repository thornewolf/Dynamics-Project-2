clear;close all;clc

% Constants
c.m = 0.25; % kg
c.k = 25; % N/m
c.l = 0.5; % m
c.I = (1/3)*c.m*c.l^2; %kg.m^2
c.g = 9.81; % m/s^2

syms l g m k I theta phi thetadot thetaddot phidot phiddot ...
    An At Cn Ct

% Cartesian locations of bar ends 
b = l*[sin(theta);-cos(theta)];
d = l*[sin(phi)+1;-cos(phi)];
rbd = d - b;
rbdSquared = rbd.^2;
% Spring Length (mag) 
ls = simplify(rbdSquared(1) + rbdSquared(2))^(1/2);
% Spring Direction (Left and Right)
ehatsL = (d - b) / ls;
ehatsR = (b - d) / ls;
% Normal / Tangential for Theta
ehatnT = [-sin(theta);cos(theta)];
ehattT = [cos(theta);sin(theta)];
% Normal / Tangential for Phi
ehatnP = [-sin(phi);cos(phi)];
ehattP = [cos(phi);sin(phi)];
% Spring Force (Left and Right)
FsL = k*(ls-l)*ehatsL;
FsR = k*(ls-l)*ehatsR;

% Spring force in normal direction (Theta)
Fs_nT = simplify(FsL(1)*ehatnT(1)+FsL(2)*ehatnT(2));
% Spring force in tangential direction (Theta)
Fs_tT = simplify(FsL(1)*ehattT(1)+FsL(2)*ehattT(2));

% Spring force in normal direction (Phi)
Fs_nP = simplify(FsR(1)*ehatnP(1)+FsR(2)*ehatnP(2));
% Spring force in tangential direction (Phi)
Fs_tP = simplify(FsR(1)*ehattP(1)+FsR(2)*ehattP(2));

% Sum of Normal Forces, left bar (Theta)
eqn(1) = m*(l/2)*thetadot^2 == Fs_nT + An - m*g*cos(theta);
% Sum of Tangential Forces, left bar (Theta)
eqn(2) = m*(l/2)*thetaddot == Fs_tT + At - m*g*sin(theta);
% Sum of Moments, left bar (Theta)
eqn(3) = -(l/2)*(m*g*sin(theta)) + (Fs_tT)*(l) == (1/3)*m*(l^2)*thetaddot;

% Sum of Normal Forces, right bar (Phi)
eqn(4) = m*(l/2)*phidot^2 == Fs_nP + Cn - m*g*cos(phi);
% Sum of Tangential Forces, right bar (Phi)
eqn(5) = m*(l/2)*phiddot == Fs_tP + Ct - m*g*sin(phi);
% Sum of Moments, right bar (Phi)
eqn(6) = -(l/2)*(m*g*sin(phi)) + (Fs_tP)*(l) == (1/3)*m*(l^2)*phiddot;

x = solve(eqn,[thetaddot,phiddot,An,At,Cn,Ct]);

%----------------------- Linearization Equations -------------------------%
LinEOMs = subs([x.thetaddot,x.phiddot],...
    {str2sym('sin(theta)'),str2sym('cos(theta)'),str2sym('sin(phi)'),str2sym('cos(phi)') ...
    str2sym('sin(phi - theta)'),str2sym('thetadot^2'),str2sym('phidot^2')},...
    {'theta',1,'phi',1,str2sym('phi - theta'),0,0});

LinEOMs = simplify(LinEOMs,'IgnoreAnalyticConstraints',true);

% Stiffness Matrix (Coefficients of Theta and Phi)
a = equationsToMatrix(LinEOMs,[theta phi]);
% Substitute Constants
a = subs(a,{'k','g','l','m'},{c.k,c.g,c.l,c.m});
% Mass Matrix
mass = [c.m 0;0 c.m];
% 
[V,D] = eig(mass\a);
% 
initialTP = [pi/36;pi/12];
syms time c1 c2 c3 c4
eqnL1 = initialTP == c1*V(:,1)*(cos(time*sqrt(D(1,1))) + c2*sin(time*sqrt(D(1,1)))) +... 
        V(:,2)*(c3*cos(time*sqrt(D(2,2))) + c4*sin(time*sqrt(D(2,2))));
eqnL2 = diff(eqnL1,time) == [0;0];

eqnL1 = subs(eqnL1, time, 0);
eqnL2 = subs(eqnL2, time, 0);
consts = solve([eqnL1,eqnL2],[c1,c2,c3,c4]);
%-------------------------------------------------------------------------%

syms theta(t) thetadot(t) phi(t) phidot(t)

thetaEOM = subs(x.thetaddot, {'theta', 'thetadot'}, {theta, thetadot});
phiEOM = subs(x.phiddot, {'phi', 'phidot'}, {phi, phidot});

LinEOMs = subs(LinEOMs, {'theta', 'thetadot','phi', 'phidot'},...
    {theta, thetadot,phi, phidot});

eom = odeFunction([thetadot; thetaEOM; phidot; phiEOM],...
    [theta; thetadot; phi; phidot],m,k,l,I,g);

eomL = odeFunction([thetadot; LinEOMs(1); phidot; LinEOMs(2)],...
    [theta; thetadot; phi; phidot],m,k,l,I,g);

theta_o = [pi/12 -pi/12 pi/36];

for i = 1:3
    [T,S] = ode45(@(t,s)eom(t,s,c.m,c.k,c.l,c.I,c.g),linspace(0,10,1001),...
        [theta_o(i),0,pi/12,0]);
    [TL,SL] = ode45(@(t,s)eomL(t,s,c.m,c.k,c.l,c.I,c.g),linspace(0,10,1001),...
        [theta_o(i),0,pi/12,0]);
    figure
    hold on
    grid on
    title('Non-Linearized')
    xlabel('Time, sec')
    ylabel('Angular Position, rad')
    plot(T,S(:,1),'DisplayName', ['\theta_o = ' num2str(rad2deg(theta_o(i))) '^o'])
    plot(T,S(:,3),'DisplayName', '\phi_o = 15^o')
    plot(TL,SL(:,1),'DisplayName', ['Linear \theta_o = ' num2str(rad2deg(theta_o(i))) '^o'])
    plot(TL,SL(:,3),'DisplayName', 'Linear \phi_o = 15^o')
    legend('show')
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [1 1 6 2.5]);
    fig = gcf;
    print(['BestFitFigure' num2str(i)],'-dpdf');
end
% 
% 
% % Checking EOM Units
% u = symunit;
% m = m*u.kg;
% g = g*u.m/u.s^2;
% l = l*u.m;
% k = k*u.N/u.m;
% An = An*u.N;
% At = At*u.N;
% Cn = Cn*u.N;
% Ct = Ct*u.N;
% theta = 'theta';
% thetadot = 'thetadot'/u.s;
% thetaddot = 'thetaddot'/u.s^2;
% phi = 'phi';
% phidot = 'phidot'/u.s;
% phiddot = 'phiddot'/u.s^2;
% 
% eqn = subs(eqn);
% unitCheck = checkUnits(eqn)

