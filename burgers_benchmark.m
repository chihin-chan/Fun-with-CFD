%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% Numerical schemes comparison with various initial conditions
% 
% Conservation law: Burger's equation
%
% By Chan Chi Hin, 12/12/19
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all 
close all

dx = 0.05;  %fixed cfl at 0.5 and decrease dx to improve solution
cfl = 0.1;  %sweet spot of CFL = 0.5, lower diffusive, higher dispersive (.6 with oscillations)
dt = cfl*dx;

%1 - Gaussian
%2 - N Function
%3 - Inverted N Function
%4 - Staircase
initial_condition = input('Select Initial Condition? \n 1. Gaussian \n 2. N Function \n 3. Inverted N \n 4. Staircase \n');

%1 - Lax-FriedRichs
%2 - Lax Wendroff
%3 - Richtmyer
%4 - MacCormack
scheme = input('Numerical Scheme?\n 1. Lax-Friedrichs \n 2. Lax-Wendroff \n 3. Richtmyer \n 4. MacCormack\n');

x=-3:dx:3;
flux = @(u) 0.5*u.^2;

%% Switching Initial Condition

u = 0*ones(length(x),1);
u_old = 0*ones(length(x),1);
u_new = 0*ones(length(x),1);
u_mid = 0*ones(length(x)-1,1);


switch initial_condition
    case 1
        for i=2:length(x)-1;
                u_old(i) = 2*exp(-x(i)^2);
                u(i) = 2*exp(-x(i)^2);
        end
    case 2
        for i=1:length(x)
            if x(i)>= -1 & x(i)<=1
                u_old(i) = 2;
                u(i) = 2;
            end
        end
    case 3
        for i=1:length(x)
            if x(i)>= -1 & x(i)<=1
                u_old(i) = 0;
                u(i) = 0;
            else 
                u_old(i) = 2;
                u(i) = 2;
            end
        end
     case 4
        for i=1:length(x)
            if x(i)>= -2 & x(i)<=-1.5
                u_old(i) = 2.5;
                u(i) = 2.5;
            elseif x(i)>-1 & x(i)<=-0.5
                u_old(i) = 2;
                u(i) = 2;
            elseif x(i)>0 & x(i)<=0.5
                u_old(i) = 1.5;
                u(i) = 1.5;
            elseif x(i)>1 & x(i)<=1.5
                u_old(i) = 1;
                u(i) = 1;
            elseif x(i)>2 & x(i)<= 2.5
                u_old(i) = 0.5;
                u(i) = 0.5;
            end
        end
end

%% Initialising Time-Step
t_final = 10;
t=0;

%% Solving
while t < t_final
    switch scheme
        %%LAX FRIEDRICHS
        case 1
            for i=2:length(x)-1
                u_new(i) = 0.5 * (u_old(i-1)+u_old(i+1)) - dt / (2*dx) * (flux(u_old(i+1))-flux(u_old(i-1)));
                num_scheme = 'Lax Friedrichs';
            end
            
        %%LAX WENDROFF    
        case 2
        for i=2:length(x)-1
            a_1 = 0.5 * ( u_old(i+1) + u_old(i) );
            a_2 = 0.5 * ( u_old(i) + u_old(i-1) );
            u_new(i) = u_old(i) - dt/(2*dx)*( flux(u_old(i+1)) - flux(u_old(i-1)) ) + ... 
                dt^2/(2*dx^2)*( a_1*( flux(u_old(i+1)) - flux(u_old(i)) ) - a_2*( flux(u_old(i)) - flux(u_old(i-1)) ));
            num_scheme = 'Lax Wendroff';
        end
        
        %Richtmyer 2-step Lax Ww
        case 3
            for i=1:length(x)-1
                u_mid(i) = 0.5*(u_old(i)+u_old(i+1))-dt/(2*dx)*(flux(u_old(i+1))-flux(u_old(i)));
                if i > 1
                    u_new(i) = u_old(i) - dt/dx*(flux(u_mid(i))-flux(u_mid(i-1)));
                end
            end
            num_scheme = 'Ritchmyer';
            
        %MacCormack
        case 4
            %Predictor
            for i=1:length(x)-1
                u_pred(i) = u_old(i) - dt/dx*(flux(u_old(i+1))-flux(u_old(i)));
            end
            %Corrector
            for i=2:length(x)-1
                u_corr(i) = u_pred(i) - dt/dx*(flux(u_pred(i))-flux(u_pred(i-1)));
            end
            %Averaging Solutions
            for i=2:length(x)-1
                u_new(i) = 0.5*(u_pred(i-1) + u_corr(i));
            end
            num_scheme = "MacCormack";
    end
    
    % Periodic Boundary Conditions
    u_new(end) = u_new(end-1); 
    u_new(1) = u_new(end);
    u_old=u_new;
    
    
    % Post-Proc
    plot(x,u, x,u_new,'-o');
    title(num_scheme);
    axis([-3 3 -1 3]);
    grid on
    drawnow
    if (t==0)
            title('PRESS <ENTER> TO START THE SIMULATION')
            pause();
    end
    t = t+dt;
end


