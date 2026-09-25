clear all
close all

%% Figure 2
% Constants
A = [1 ; 1]; % a_s, a_r
B = [0.3 0.25
    0.25 0.3]; % interaction matrix (where bji is the effect of species i on j)
x0 = [0.5 ; 1]; % initial conditions (unless specified) [x_s, x_r]
tBounds = [0, 1000]; % time interval to consider


figure(1) % Figure 2A
plotLV(tBounds, x0, B, A, {"2A:","Coexistence in Absence of a Pulse"});
figure(2) % Figure 2B
plot_pulsed_2d_bifurcation_coexistence(A, B);

%% Figure 3
% Constants
A = [1 ; 1]; % a_s, a_r
B = [0.3 0.25
    0.4 0.3]; % interaction matrix (where bji is the effect of species i on j)
x0 = [0.5 ; 0.5]; % initial conditions (unless specified) [x_s, x_r]
tBounds = [0, 1000]; % time interval to consider

figure(3) % Figure 3A
plotLV(tBounds, x0, B, A, {"3A:","Sensitive Cells Have Competitive Advantage in the Absence of a Pulse"});
figure(4) % Figure 3B
plot_pulsed_2d_bifurcation_competitive_disadvantage(A, B);

A = [1 ; 0.5]; % a_s, a_r
B = [0.3 0.25
    0.25 0.3]; % interaction matrix (where bji is the effect of species i on j)

figure(5) % not included in manuscript 
plotLV(tBounds, x0, B, A, {"3(D):","Sensitive Cells Have Growth Advantage in the Absence of a Pulse"});

h_star = (1- sqrt(5))/2; tau_star = - 1 * log(1 + h_star) / A(1);
figure(6) % Figure 3C
plot_pulsed_2d_bifurcation_growth_disadvantage(A, B, h_star, tau_star);

%% Figure 4
% Constants
A = [1 ; 1]; % a_s, a_r
B = [0.3 0.4
    0.25 0.3]; % interaction matrix (where bji is the effect of species i on j)
x0 = [0.5 ; 0.5]; % initial conditions (unless specified) [x_s, x_r]
tBounds = [0, 1000]; % time interval to consider

figure(7) % Figure 4A 
plotLV(tBounds, x0, B, A, {"4A:","Resistant Cells Have Competitive Advantage in Absence of a Pulse"});
figure(8) % Figure 4B
plot_pulsed_2d_bifurcation_competitive_advantage(A, B);

A = [0.5 ; 1]; % a_s, a_r
B = [0.3 0.25
    0.25 0.3];
figure(9) % not included in manuscript 
plotLV(tBounds, x0, B, A, {"4(C):","Resistant Cells Have Growth Advantage in Absence of a Pulse"});
figure(10) % Figure 4D
plot_pulsed_2d_bifurcation_competitive_advantage(A, B);
hold on;
title({"4(D):","Resistance Confers Competitive Advantage","in Untreated System"}, FontSize=18, Interpreter="latex");
hold off;

%% Figure 5
% Constants
A = [-1 1]; % a_x, a_y
B= [0 -1/2;
    1/2 1/3]; % interaction matrix (where bji is the effect of species i on j)

figure(11); % Figure 5
plot_pulsed_swd_bifurcation(A, B);

%% Figure 6
% example differential system
function dx = test_system_1(t, x)
    dx = x;
end

figure(12); % Figure 6A
plot_pulsed([0,2.01],3,@test_system_1,2,-0.75,{"6A:","Pulsed Nonautonomous System"});
figure(13); % Figure 6B
plot_periodic([0,3],3,@test_system_1,2,-0.75,{"6B:","Continuous Autonomous System"});

%% General Functions for Plotting.

% L.V. System Definition and Plotting
function dx = LVSystem(t, x, D, A)
dx = (A - D*x).*x;
end

% Generate labels for plotting
function labs = getLabels(n)
labs = cell(n,1);
for i = 1:n
    labs{i} = "Species " + i;
end
end

% Plot a given LV system
function [tode,yode] = plotLV(tBounds, xInit, D, A, name)

% Generate numerical solution
[tode,yode] = ode45(@(t,y) LVSystem(t,y,D,A), tBounds, xInit);

plot(tode, yode(:,1), "LineStyle","-","LineWidth",2);
hold on;
for i = 2:length(yode(1,:))
    plot(tode, yode(:,i), "LineStyle","--","LineWidth",2);
end

xlim(tBounds);
ylim([0,max(yode, [], "all") + 0.6]);
xlabel("Time", "FontSize",14,"FontName","Cambria", "Interpreter","latex");
ylabel("Population Size", "FontSize",14,"FontName","Cambria","Interpreter","latex");
title(name, "FontSize", 20, "FontName","Cambria","Interpreter","latex");    

% Generate legend
legend(getLabels(size(xInit,1)),"Location","north","NumColumns",2)
hold off;
end

% CASE ONE: COEXISTENCE IN UNPULSED SYSTEM
% resistance does not confer an advantage or disadvantage
function plot_pulsed_2d_bifurcation_coexistence(a, b)
% plot styling
plot(linspace(-1,1),zeros(1,100), LineWidth=1.25, Color="#5c5c5c");
hold on;
plot(zeros(1,100), linspace(0,10), LineWidth=1.25, Color="#5c5c5c");
xlim([-1,0])
ylim([0, 5])
title({"2B:","Untreated Coexistence"}, FontSize=18, Interpreter="latex")
xlabel("Baseline Pulse Effect [$h_0$]", Interpreter="latex", FontSize=18);
ylabel("Pulse Period [$\tau$]", Interpreter="latex", FontSize=18);
hAxes = findobj(gcf,"Type","axes");
hAxes.FontSize = 18;

% parameter regions of interest    
tau = linspace(0,10);
h0 = linspace(-0.999,0);
h1 = h0;              % sensitive species
h2 = -1*h0(1:100).^2; % resistant species

K1 = (a(1) - a(2)*b(1,2)/b(2,2))^(-1); 
K2 = (a(2) - a(1)*b(2,1)/b(1,1))^(-1);
p1 = (log(1+h2).*b(1,2)/b(2,2) - log(1+h0)) * K1; % boundary of species 1 persistence  
p2 = (log(1+h0)*b(2,1)/b(1,1) - log(1+h2)) * K2; % boundary of species 2 persistence  

% priority effect - appears in tau < 0 region
patch([], [], [0.9412, 0.7725, 0.4431], "FaceAlpha", 0.54)

% coexistence
% patch([h0(1:100),0,h0(100:-1:1)],[p1(1:100),0,ones(1,100)*p1(1)],[0.6471 0.3490 0.6667], "FaceAlpha", 0.54)
patch([h0(1:100),0,h0(100:-1:1)],[p1(1:100),0,ones(1,100)*p1(1)],[1 1 1], "FaceAlpha", 0.54)

% system extinction
patch([h0(1:100),0,h0(100:-1:1)], [-1*log(h2 + 1)/a(2), 0, zeros(1,100)],[0.8078, 0.8078, 0.8078],"FaceAlpha",0.54);

% species 1 wins - does not appear in h0 < 0 region
patch([], [], [0.551 0.824 0.777], "FaceAlpha", 0.54);

% species 2 wins - where 2 is resistant to treatment
patch([h0(1:100),0,h0(100:-1:1)],[-1*log(h2 + 1)/a(2), 0, p1(100:-1:1)], [0.8784, 0.1686, 0.2078], "FaceAlpha", 0.54);

hold off;
end


% CASE TWO A: COMPETITIVE ADVANTAGE OF SENSITIVE SPECIES IN UNTREATED
% SYSTEM
% resistance confers a competitive disadvantage 
function plot_pulsed_2d_bifurcation_competitive_disadvantage(a, b)

% plot styling
plot(linspace(-1,1),zeros(1,100), LineWidth=1.25, Color="#5c5c5c");
hold on;
plot(zeros(1,100), linspace(0,10), LineWidth=1.25, Color="#5c5c5c");
xlim([-1,0])
ylim([0, 5])
title({"3B:","Resistance Confers Competitive Disadvantage","in Untreated System"}, FontSize=18, Interpreter="latex")
xlabel("Baseline Pulse Effect [$h_0$]", Interpreter="latex", FontSize=18);
ylabel("Pulse Period [$\tau$]", Interpreter="latex", FontSize=18);
hAxes = findobj(gcf,"Type","axes");
hAxes.FontSize = 18;

% parameter regions of interest    
tau = linspace(0,11);
h0 = linspace(-0.999,0);
h1 = h0;              % sensitive species
h2 = -1*h0(1:100).^2; % resistant species

K1 = (a(1) - a(2)*b(1,2)/b(2,2))^(-1); 
K2 = (a(2) - a(1)*b(2,1)/b(1,1))^(-1);
p1 = (log(1+h2).*b(1,2)/b(2,2) - log(1+h0)) * K1; % boundary of species 1 persistence  
p2 = (log(1+h0)*b(2,1)/b(1,1) - log(1+h2)) * K2; % boundary of species 2 persistence  

% priority effect - appears in tau < 0 region
patch([h0(1:100),0,h0(100:-1:1)], [p2(1:100),0,p1(100:-1:1)], [0.9412, 0.7725, 0.4431], "FaceAlpha", 0.54)

% coexistence - does not appear in tau > 0; h < 0 region
% patch([ ],[ ],[0.6471 0.3490 0.6667], "FaceAlpha", 0.54)

% system extinction
patch([h0(1:100),0,h0(100:-1:1)], [-1*log(h2 + 1)/a(2), 0, zeros(1,100)],[0.8078, 0.8078, 0.8078],"FaceAlpha",0.54);

% species 1 wins - does not appear in h0 < 0 region
patch([h0(1:100),0,h0(100:-1:1)], [p1(1:100),0,11*ones(1,100)], [0.551 0.824 0.777], "FaceAlpha", 0.54);

% species 2 wins - where 2 is resistant to treatment
patch([h0(1:100),0,h0(100:-1:1)],[-1*log(h2 + 1)/a(2), 0, p2(100:-1:1)], [0.8784, 0.1686, 0.2078], "FaceAlpha", 0.54);

hold off;
end

% CASE TWO B: GROWTH DISADVANTAGE OF RESISTANT SPECIES IN UNTREATED
% SYSTEM
function plot_pulsed_2d_bifurcation_growth_disadvantage(a, b, h_star, tau_star)
% plot styling
plot(linspace(-1,1),zeros(1,100), LineWidth=1.25, Color="#5c5c5c");
hold on;
plot(zeros(1,100), linspace(0,10), LineWidth=1.25, Color="#5c5c5c");
xlim([-1,0]);
ylim([0, 1]);
title({"3C:","Resistance Confers Growth Disadvantage","in Untreated System"}, FontSize=18, Interpreter="latex")
xlabel("Baseline Pulse Effect [$h_0$]", Interpreter="latex", FontSize=18);
ylabel("Pulse Period [$\tau$]", Interpreter="latex", FontSize=18);
hAxes = findobj(gcf,"Type","axes");
hAxes.FontSize = 18;

% parameter regions of interest    
tau = [linspace(0,tau_star) linspace(tau_star,11)];
h0 = [linspace(-0.999,h_star) linspace(h_star,0)];
h1 = h0;              % sensitive species
h2 = -1*h0.^2; % resistant species

K1 = (a(1) - a(2)*b(1,2)/b(2,2))^(-1); 
K2 = (a(2) - a(1)*b(2,1)/b(1,1))^(-1);
p1 = (log(1+h2).*b(1,2)/b(2,2) - log(1+h0)) * K1; % boundary of species 1 persistence  
p2 = (log(1+h0)*b(2,1)/b(1,1) - log(1+h2)) * K2; % boundary of species 2 persistence  

% priority effect - does not appear in tau > 0; h < 0 region
% patch([], [], [0.9412, 0.7725, 0.4431], "FaceAlpha", 0.54)

% coexistence 
% patch([h0(101:200),0,h0(200:-1:101)],[p2(101:200),0,p1(200:-1:101)],[0.6471 0.3490 0.6667], "FaceAlpha", 0.54)
patch([h0(101:200),0,h0(200:-1:101)],[p2(101:200),0,p1(200:-1:101)],[1 1 1], "FaceAlpha", 0.54)

% system extinction
patch([h0(1:100),h0(101:200),h0(200:-1:1)], [-1*log(h0(1:100)+1)/a(1),-1*log(h2(101:200) + 1)/a(2),zeros(1,200)],[0.8078, 0.8078, 0.8078],"FaceAlpha",0.54);

% species 1 wins 
patch([h0(1:100),h0(101:200),0,h0(200:-1:1)], [-1*log(h0(1:100)+1)/a(1),p2(101:200),0,11*ones(1,200)], [0.551 0.824 0.777], "FaceAlpha", 0.54);

% species 2 wins - where 2 is resistant to treatment
patch([h0(101:200),0,h0(200:-1:101)],[p1(101:200),0,-1*log(h2(200:-1:101) + 1)/a(2)], [0.8784, 0.1686, 0.2078], "FaceAlpha", 0.54);

hold off;
end

% CASE THREE A: COMPETITIVE ADVANTAGE OF RESISTANT SPECIES IN UNTREATED
% SYSTEM
function plot_pulsed_2d_bifurcation_competitive_advantage(a, b)
% plot styling
plot(linspace(-1,1),zeros(1,100), LineWidth=1.25, Color="#5c5c5c");
hold on;
plot(zeros(1,100), linspace(0,10), LineWidth=1.25, Color="#5c5c5c");
xlim([-1,0])
ylim([0, 5])
title({"4B:","Resistance Confers Competitive Advantage","in Untreated System"}, FontSize=18, Interpreter="latex")
xlabel("Baseline Pulse Effect [$h_0$]", Interpreter="latex", FontSize=18);
ylabel("Pulse Period [$\tau$]", Interpreter="latex", FontSize=18);
hAxes = findobj(gcf,"Type","axes");
hAxes.FontSize = 18;

% parameter regions of interest    
tau = linspace(0,11);
h0 = linspace(-0.999,0);
h1 = h0;              % sensitive species
h2 = -1*h0(1:100).^2; % resistant species

K1 = (a(1) - a(2)*b(1,2)/b(2,2))^(-1); 
K2 = (a(2) - a(1)*b(2,1)/b(1,1))^(-1);
p1 = (log(1+h2).*b(1,2)/b(2,2) - log(1+h0)) * K1; % boundary of species 1 persistence  
p2 = (log(1+h0)*b(2,1)/b(1,1) - log(1+h2)) * K2; % boundary of species 2 persistence  

% priority effect - appears in tau < 0 region
% patch([], [], [0.9412, 0.7725, 0.4431], "FaceAlpha", 0.54)

% coexistence - does not appear in tau > 0; h < 0 region
% patch([],[],[0.6471 0.3490 0.6667], "FaceAlpha", 0.54)

% system extinction
patch([h0(1:100),0,h0(100:-1:1)], [-1*log(h2 + 1)/a(2), 0, zeros(1,100)],[0.8078, 0.8078, 0.8078],"FaceAlpha",0.54);

% species 1 wins - does not appear in h0 < 0 region
% patch([], [], [0.551 0.824 0.777], "FaceAlpha", 0.54);

% species 2 wins - where 2 is resistant to treatment
patch([h0(1:100),0,h0(100:-1:1)],[-1*log(h2 + 1)/a(2), 0, 11*ones(1,100)], [0.8784, 0.1686, 0.2078], "FaceAlpha", 0.54);

hold off;
end

% SWD BIFURCATION DIAGRAM: EFFECT OF VARYING PULSE PERIOD V. TREATMENT SPECIFICITY
function plot_pulsed_swd_bifurcation(a, b)
% plot styling
plot(linspace(-1,1),zeros(1,100), LineWidth=1.25, Color="#5c5c5c");
hold on;
plot(zeros(1,100), linspace(0,10), LineWidth=1.25, Color="#5c5c5c");
xlim([-1,0])
ylim([0, 5])
title({"5:","SWD: Tradeoffs with Targeted Treatment"}, FontSize=18, Interpreter="latex")
xlabel("Pulse Effect on SWD [$h_d$]", Interpreter="latex", FontSize=18);
ylabel("Pulse Period [$\tau$]", Interpreter="latex", FontSize=18);
hAxes = findobj(gcf,"Type","axes");
hAxes.FontSize = 18;

% parameter regions of interest    
tau = [linspace(0,7)];
h0 = [linspace(-0.999,0,200)];
hw = -1-h0;              % predator species (1)
hd = h0; % prey species (2)

perm = -1 / (b(2,1)/b(2,2)*a(2) - abs(a(1))) * log(-hd.*(1+hd).^(b(2,1)/b(2,2))); % boundary of system permanence
pd = -1*(log(1+hd))./a(2); % boundary of SWD survival

% priority effect - does not appear in tau > 0; h < 0 region
% patch([], [], [0.9412, 0.7725, 0.4431], "FaceAlpha", 0.54)

% permanence 
% patch([h0(1:199),h0(200),h0(200:-1:1)],[perm(1:199),perm(1),perm(1)*ones(1,200)],[0.6471 0.3490 0.6667], "FaceAlpha", 0.54)
patch([h0(1:199),h0(200),h0(200:-1:1)],[perm(1:199),perm(1),perm(1)*ones(1,200)],[1 1 1], "FaceAlpha", 0.54)

% system extinction
patch([h0(1:200),h0(200:-1:1)], [pd(1:200),zeros(1,200)],[0.8078, 0.8078, 0.8078],"FaceAlpha",0.54);

% SWD 
patch([h0(1:200),h0(199:-1:1)], [pd(1:200),perm(199:-1:1)], [0.551 0.824 0.777], "FaceAlpha", 0.54);

% wasp wins 
% patch([],[], [0.8784, 0.1686, 0.2078], "FaceAlpha", 0.54);

end

% NUMERICAL SOLUTION TO PIECEWISE CONTINUOUS SYSTEM
function [tode,yode] = plot_periodic(tBounds, xInit, f, tau, h, name)
% Generate numerical solution
tInit = tBounds(1); 
subInt = [tInit,min(tBounds(1)+tau, tBounds(2))]; % current subinterval being considered
[tode, yode] = ode45(@(t,y) f(t,y), subInt, xInit); % normal growth solution

k = floor(tBounds(2)/tau);
while (subInt(2) < tBounds(2)) 
    % apply the pulse to the last time step
    tInit = tode(length(tode)); xInit = yode(length(yode),:); % initial condition before the pulse
    subInt = [tInit, min(tInit+1, tBounds(2))]; % new subinterval
    [sub_tode, sub_yode] = ode45(@(t,y) y.*log(1+h), subInt, xInit); % solution

    tode = [tode; sub_tode]; yode = [yode; sub_yode];

    % check if the while condition is still met
    if(subInt(2) >= tBounds(2)) break; end

    % % advance to the next subinterval with the new initial condition
    tInit = tode(length(tode));  xInit = yode(length(yode),:); % new initial condition
    subInt = [tInit, min(tInit+tau, tBounds(2))]; % new subinterval

    % calculate the solution and add it to the solution vector
    [sub_tode, sub_yode] = ode45(@(t,y) f(t,y), subInt, xInit); % solution
    tode = [tode; sub_tode]; yode = [yode; sub_yode];
end

plot(tode, yode(:,1), "LineStyle","-","LineWidth",2);
hold on;
for i = 2:length(yode(1,:))
    plot(tode, yode(:,i), "LineStyle","--","LineWidth",2);
end

xlim(tBounds);
ylim([0,max(yode, [], "all") + 0.6]);
xlabel("Time", "FontSize",14,"FontName","Cambria");
ylabel("Population Size", "FontSize",14,"FontName","Cambria");
title(name, "FontSize", 20, "FontName","Cambria");    

% Generate legend
legend(getLabels(length(yode(1,:))),"Location","northwest","NumColumns",2)
hold off;
end

% NUMERICAL SOLUTION TO A PULSED SYSTEM
% Plot a pulsed L.V. system
function [tode,yode] = plot_pulsed(tBounds, xInit, f, tau, h, name)
% Generate numerical solution
tInit = tBounds(1); 
subInt = [tInit,min(tBounds(1)+tau, tBounds(2))]; % current subinterval being considered
[tode, yode] = ode45(@(t,y) f(t,y), subInt, xInit); % solution

while (subInt(2) < tBounds(2)) 
    % apply the pulse to the last time step
    yode = [yode; yode(length(yode),:) .* (1 + h)];
    tode = [tode; tode(length(tode))];

    % % advance to the next subinterval with the new initial condition
    tInit = tode(length(tode));  xInit = yode(length(yode),:); % new initial condition
    subInt = [tInit, min(tInit+tau, tBounds(2))]; % new subinterval

    % calculate the solution and add it to the solution vector
    [sub_tode, sub_yode] = ode45(@(t,y) f(t,y), subInt, xInit); % solution
    tode = [tode; sub_tode]; yode = [yode; sub_yode];
end

plot(tode, yode(:,1), "LineStyle","-","LineWidth",2);
hold on;
for i = 2:length(yode(1,:))
    plot(tode, yode(:,i), "LineStyle","--","LineWidth",2);
end

xticks(tBounds(1):tau:tBounds(2));
xticklabels(tBounds(1):tau:tBounds(2));
xlim(tBounds);
ylim([0,max(yode, [], "all") + 0.6]);
xlabel("Time", "FontSize",14,"FontName","Cambria");
ylabel("Population Size", "FontSize",14,"FontName","Cambria");
title(name, "FontSize", 20, "FontName","Cambria");    

% Generate legend
legend(getLabels(length(yode(1,:))),"Location","northwest","NumColumns",2);
hold off;
end

