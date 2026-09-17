clear all
close all

%% Figure 2
% Constants
A = [1 ; 1]; % a_s, a_r
B = [0.3 0.25
    0.25 0.3]; % interaction matrix (where bji is the effect of species i on j)
x0 = [0.5 ; 1]; % initial conditions (unless specified) [x_s, x_r]
tBounds = [0, 1000]; % time interval to consider

figure(1)
plotLV(tBounds, x0, B, A, {"Case A:","permanent in the absence of a pulse"});
figure(2)
plot_pulsed_2d_bifurcation_coexistence(A, B);

%% Figure 3
A = [1 ; 1]; % a_s, a_r
B = [0.3 0.25
    0.4 0.3]; % interaction matrix (where bji is the effect of species i on j)
x0 = [0.5 ; 0.5]; % initial conditions (unless specified)
tBounds = [0, 1000]; % time interval to consider

figure(3)
plotLV(tBounds, x0, B, A, {"Case A:","permanent in the absence of a pulse"});
figure(4)
plot_pulsed_2d_bifurcation_competitive_disadvantage(A, B)
ylim([0, 5])

A = [1 ; 0.5]; % a_s, a_r 
B = [0.3 0.25
    0.25 0.3]; % interaction matrix (where bji is the effect of species i on j)
figure(5)
h_star = (1- sqrt(5))/2; tau_star = - 1 * log(1 + h_star) / A(1);
plot_pulsed_2d_bifurcation_growth_disadvantage(A, B,h_star, tau_star);
ylim([0,1])

% General L.V. System Definition and Plotting
% Define the differential system
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

% Plot a given system
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

% Figure 2: A
% resistance does not confer an advantage or disadvantage
function plot_pulsed_2d_bifurcation_coexistence(a, b)
    % plot styling
    plot(linspace(-1,1),zeros(1,100), LineWidth=1.25, Color="#5c5c5c");
    hold on;
    plot(zeros(1,100), linspace(0,10), LineWidth=1.25, Color="#5c5c5c");
    xlim([-1,0])
    ylim([0, 5])
    title("Pulsed Treatment Bifurcation Plot", FontSize=18, Interpreter="latex")
    xlabel("Baseline Pulse Effect [$$h_0$$]", Interpreter="latex", FontSize=16);
    ylabel("Pulse Period [$$\tau$$]", Interpreter="latex", FontSize=16);
    
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
    patch([], [], [.68 .6 .01], "FaceAlpha", 0.54)
    
    % coexistence
    patch([h0(1:100),0,h0(100:-1:1)],[p1(1:100),0,ones(1,100)*p1(1)],[.56 .69 .95], "FaceAlpha", 0.54)
    
    % system extinction
    patch([h0(1:100),0,h0(100:-1:1)], [-1*log(h2 + 1)/a(2), 0, zeros(1,100)],[0.5, 0.691,0.824],"FaceAlpha",0.54);
    
    % species 1 wins - does not appear in h0 < 0 region
    patch([], [], [0.551 0.824 0.777], "FaceAlpha", 0.54);
    
    % species 2 wins - where 2 is resistant to treatment
    patch([h0(1:100),0,h0(100:-1:1)],[-1*log(h2 + 1)/a(2), 0, p1(100:-1:1)], [.98 .5 0.445], "FaceAlpha", 0.54);

end

%% Figure 3B
function plot_pulsed_2d_bifurcation_competitive_disadvantage(a, b)
    % plot styling
    plot(linspace(-1,1),zeros(1,100), LineWidth=1.25, Color="#5c5c5c");
    hold on;
    plot(zeros(1,100), linspace(0,10), LineWidth=1.25, Color="#5c5c5c");
    xlim([-1,0])
    ylim([0, 10])
    title("Resistance Confers Competitive Disadvantage" + newline+ "in Untreated System", FontSize=18, Interpreter="latex")
    xlabel("Baseline Pulse Effect [$$h_0$$]", Interpreter="latex", FontSize=16);
    ylabel("Pulse Period [$$\tau$$]", Interpreter="latex", FontSize=16);
    
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
    patch([h0(1:100),0,h0(100:-1:1)], [p2(1:100),0,p1(100:-1:1)], [.68 .6 .01], "FaceAlpha", 0.54)
    
    % coexistence - does not appear in tau > 0; h < 0 region
    % patch([h0(1:100),0,h0(100:-1:1)],[p1(1:100),0,ones(1,100)*p1(1)],[.56 .69 .95], "FaceAlpha", 0.54)
    
    % system extinction
    patch([h0(1:100),0,h0(100:-1:1)], [-1*log(h2 + 1)/a(2), 0, zeros(1,100)],[0.5, 0.691,0.824],"FaceAlpha",0.54);
    
    % species 1 wins - does not appear in h0 < 0 region
    patch([h0(1:100),0,h0(100:-1:1)], [p1(1:100),0,11*ones(1,100)], [0.551 0.824 0.777], "FaceAlpha", 0.54);
    
    % species 2 wins - where 2 is resistant to treatment
    patch([h0(1:100),0,h0(100:-1:1)],[-1*log(h2 + 1)/a(2), 0, p2(100:-1:1)], [.98 .5 0.445], "FaceAlpha", 0.54);
end

% Figure 3 C: GROWTH DISADVANTAGE OF RESISTANT SPECIES IN UNTREATED
% SYSTEM
% resistance confers a competitive disadvantage 
function plot_pulsed_2d_bifurcation_growth_disadvantage(a, b, h_star, tau_star)
    % plot styling
    plot(linspace(-1,1),zeros(1,100), LineWidth=1.25, Color="#5c5c5c");
    hold on;
    plot(zeros(1,100), linspace(0,10), LineWidth=1.25, Color="#5c5c5c");
    xlim([-1,0])
    ylim([0, 10])
    title("Pulsed Treatment Bifurcation Plot", FontSize=18, Interpreter="latex")
    xlabel("Baseline Pulse Effect [$$h_0$$]", Interpreter="latex", FontSize=16);
    ylabel("Pulse Period [$$\tau$$]", Interpreter="latex", FontSize=16);
    
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
    % patch([], [], [.68 .6 .01], "FaceAlpha", 0.54)
    
    % coexistence 
    patch([h0(101:200),0,h0(200:-1:101)],[p2(101:200),0,p1(200:-1:101)],[.56 .69 .95], "FaceAlpha", 0.54)
    
    % system extinction
    patch([h0(1:100),h0(101:200),h0(200:-1:1)], [-1*log(h0(1:100)+1)/a(1),-1*log(h2(101:200) + 1)/a(2),zeros(1,200)],[0.5, 0.691,0.824],"FaceAlpha",0.54);
    
    % species 1 wins 
    patch([h0(1:100),h0(101:200),0,h0(200:-1:1)], [-1*log(h0(1:100)+1)/a(1),p2(101:200),0,11*ones(1,200)], [0.551 0.824 0.777], "FaceAlpha", 0.54);
    
    % species 2 wins - where 2 is resistant to treatment
    patch([h0(101:200),0,h0(200:-1:101)],[p1(101:200),0,-1*log(h2(200:-1:101) + 1)/a(2)], [.98 .5 0.445], "FaceAlpha", 0.54);
end



