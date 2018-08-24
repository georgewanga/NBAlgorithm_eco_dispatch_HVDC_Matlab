clc; % clear any work or data in the command window
clear all; % clear all varriable values before use
close all; % close all open figures

doc_name = 'ED_result.doc';
plot_Fcost = 'FuelCostCurve.png';
plot_Iterr = 'ItterationsCurve.png';
plot_Ploss = 'PowerLossCurve.png';
bar_Ploss = 'PowerLossChart.png';
bar_Fcost = 'FuelCostBar.png';
hvdc_Losses = 'HVDC_loses.png';
transmission_modes = ["HVAC","HVDC"];
source = ["6thermal","4thermal","2wind"];

% prelocating matrices that change in length
[power_loss,F_cost,iterrations,sw_loss,cond_loss,tl_loss,F_cost_inst] = deal(zeros);

demand = [120 150 180 210 240 270 300 330 360 390 420]; % load demands
load_demand_values = numel(demand); % numel counts the elements of matrix

print = fopen(doc_name,'w+');

% variables available to all functions
global fuel_coefficients B power_demand Pg_limits transmission_type ...
    Cond_loss SW_loss TL_loss convergence_time start_timing DRi URi ...
    n f_cost beta tao time instability inst_const

% fuel_coefficients matrix having 5 columns of fuel cost coefficients
fuel_coefficients = [0.00375 2.00 240  0 0;
    0.01750 1.75 200 0 0;
    0.06250 1.00 220 40 0.008;
    0.00834 3.25 200 30 0.009;
    0.02500 3.00 220 0 0;
    0.02500 3.00 190 0 0];
generator_limits = [50 200;20 80;15 50;10 35;10 30;12 40];
%RAMP RATE CONSTRAINTS
DRi= [85 22 15 16 9 16];
URi= [65 12 12 8 6 8];
beta = 1.75;
tao = 2.85;
time = 10; % instability time in seconds
instability = false(); % set the first calculations to be without instability
n = length(fuel_coefficients(:,1)); %Returns the length of the fuel_coefficients variable

for type = 1:numel(transmission_modes)% looping through each mode
    transmission_type = transmission_modes(type);
    fprintf(print,strcat('ECONOMIC DISPATCH FOR _',transmission_type,...
        ' USING NOVEL BAT OPTIMIZATION ALGORITHM \n'));
    %% Step 1:finding the B matrix
    loss_coef = [0.000218 0.000103 0.000009 -0.000010 0.000002 0.000027
        0.000103 0.000181 0.000004 -0.000015 0.000002 0.000030
        0.000009 0.000004 0.000417 -0.000131 -0.000153 -0.000107
        -0.000010 -0.000015 -0.000131 0.000221 0.000094 0.000050
        0.000002 0.000002 -0.000153 0.000094 0.000243 -0.000000
        0.000027 0.000030 -0.000107 0.000050 -0.000000 0.000358];
    %% Step 2: getting power demand and setting incremental cost(lamda)
    for idx = 1:load_demand_values
        power_demand = demand(idx);
        disp(strcat('Computing dispatch for >',num2str(power_demand),...
            'MW in >',transmission_type,', ',num2str(load_demand_values...
            -idx),' more values to go...'))
        disp('Working please wait ...')
        %% Step3: Deploying Novel Bat Algorithm (NBA)
        if (min(generator_limits(:,1)) <= power_demand)&&(power_demand <= sum(generator_limits(:,2)))
            % setting the parameters in the basic Novel Bat Algorithm (NBA)
            M = 1000;   %number of iterations
            pop = 30;
            gamma = 0.9;
            alpha = 0.99;
            r0Max = 1;
            r0Min = 0;
            AMax = 2;
            AMin = 1;
            freqDMax = 1.5;
            freqDMin = 0;

            % setting the additional parameters in Novel Bat Algorithm (NBA)
            G = 10;
            probMax = 0.9;
            probMin = 0.6;
            thetaMax = 1;
            thetaMin = 0.5;
            wMax = 0.9;
            wMin = 0.5;
            CMax = 0.9;
            CMin = 0.1;
            if strcmp(transmission_type,'HVDC') %assigning a different B for HVDC
                B = 0.45*loss_coef;
            else
                B = loss_coef;
            end
            % starting a timer and saving the time to an output argument start_timing
            start_timing = tic;
            % calling the Novel Bat Optimization algorithm
            dim = length(B(:,1));
            Pg_limits = generator_limits;
            [bestX, bestIter] = NBA( @eldnba, M, pop, dim, G, gamma, ...
                alpha, r0Max, r0Min, AMax, AMin, freqDMax, freqDMin, probMax, ...
                probMin, CMax, CMin, thetaMax, thetaMin, wMax, wMin );
            % toc work together with tic to measure elapsed time
            convergence_time = toc(start_timing);
            %% Step 4: calculations and results display
            [fuel_cost, gen_powers, P_loss] = eldnba(bestX);
            fuel_cost_inst = fuel_cost + inst_const; 
            error = sum(gen_powers) - power_demand - P_loss;

            fprintf(print,'\nPower Demand = %g MW',power_demand);

            for gens = 1:n
                fprintf(print,'\nPower from Generator %i = %g MW',gens,gen_powers(gens));
            end
            fprintf(print,'\n Total Power Generated = %g MW',sum(gen_powers));
            if strcmp(transmission_type,'HVDC') 
                fprintf(print,'\nConduction Power Loss = %g MW',Cond_loss);
                fprintf(print,'\nSwitching Power Loss = %g MW',SW_loss);
            end
            fprintf(print,'\nTransmission Power Loss = %g MW',TL_loss);
            fprintf(print,'\n Total Power Loss = %g MW',P_loss);
            for gens = 1:n
                fprintf(print,'\nFuel cost for Generator %i = $%g',gens,f_cost(gens));
            end
            fprintf(print,'\n Total Fuel cost = $%g',fuel_cost);
            if strcmp(transmission_type,'HVDC')
                fprintf(print,'\n Total Fuel cost with instability = $%g',fuel_cost_inst);
            end
            fprintf(print,'\nError = %g MW', error);
            fprintf(print,'\nNumber of iterrations = %g',bestIter);
            fprintf(print,'\nConvergence time = %g seconds',convergence_time);
            fprintf(print,'\n\n********************************************\n');
            % create matrix values for plotting graphs
            if strcmp(transmission_type,'HVAC') %assigning a different B for HVDC
                power_loss(idx) = P_loss; % power loss matrix
                F_cost(idx) = fuel_cost; % total fuel cost matrix
                counter = idx;
            else
                power_loss(counter+idx) = P_loss; % power loss matrix
                F_cost(counter+idx)= fuel_cost; % total fuel cost matrix
                F_cost_inst(idx)= fuel_cost_inst; % total fuel cost matrix
                iterrations(idx) = bestIter; %number of itterations matrix
                sw_loss(idx) = SW_loss;
                cond_loss(idx) = Cond_loss;
                tl_loss(idx) = TL_loss;
            end
        else
            disp(['Edit demands for values between ',num2str(min(...
                Pg_limits(:,1))),'MW and ',num2str(sum(Pg_limits(:,2))),'MW'])
        end % closing the if analysis
    end     % jump out of the no of load demands for loop
end     % jump out of the no of transmission modes for loop
fclose(print);
disp('Ploting on Graphs, hold on ...')

fig1 = figure('Visible', 'off');
% plot HVAC
plot(demand,F_cost(1,1:load_demand_values),'--rs');
hold on; % KEEPS THE AXES ACTIVE FOR NEXT PLOT
% plot HVDC
plot(demand,F_cost(1,load_demand_values+1:load_demand_values*2),'--bo');
hold on; % KEEPS THE AXES ACTIVE FOR NEXT PLOT
% plot HVDC
plot(demand,F_cost_inst(1,1:load_demand_values),'--go');
hold off; % closes the axes
% set a title and key (legend) for the figure
title('Graph of Fuel Cost(6 thermal generators) against Power Demand'); %title
legend(transmission_modes(1),strcat(transmission_modes(2),' stable'), ...
    strcat(transmission_modes(2),' with instability'),'Location','northwest');
%label axes
xlabel('Power Demand (MW)'); % x-axis label
ylabel('Fuel Cost ($)'); % y-axis label
grid on;                    %enable grid on the graph
saveas(fig1,plot_Fcost);    % saves the graph to file

fig2 = figure('Visible', 'off');
plot(F_cost(1,load_demand_values+1:load_demand_values*2),iterrations,'--rs','LineWidth',2,...
    'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',7);
% set a title and key (legend) for the figure
title('Fuel Cost against Iterrations(6 thermal gen) for ED with HVDC'); %title
legend('HVDC','Location','northeast');
%label axes
xlabel('Fuel Cost ($)'); % x-axis label
ylabel('Number of Iterrations'); % y-axis label
grid on;                    %enable grid on the graph
saveas(fig2,plot_Iterr);    % saves the graph to file

fig3 = figure('Visible', 'off');
y4 = [power_loss(1,1:load_demand_values);power_loss(1,load_demand_values+1:load_demand_values*2)]';
bar(demand,y4);
% set a title and key (legend) for the figure
title('Graph of Power Loss (6 thermal gen) against Power Demand');
%label axes
xlabel('Power Demand (MW)'); % x-axis label
ylabel('Power Loss (MW)');  % y-axis label
grid on;                    %enable grid on the graph
legend(transmission_modes(1),transmission_modes(2),'Location','northwest');
%legend(p(:,:),{'HVDC','HVAC'},'Location','north');
saveas(fig3,bar_Ploss);    % saves the graph to file

fig4 = figure('Visible', 'off');
% plot HVAC
plot(demand,power_loss(1,1:load_demand_values),'--rs','LineWidth',2,...
    'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',7);
hold on; % KEEPS THE AXES ACTIVE FOR NEXT PLOT
% plot HVDC
plot(demand,power_loss(1,load_demand_values+1:load_demand_values*2),'--bo','LineWidth',2,...
    'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',7);
hold off; % closes the axes
legend(transmission_modes(1),transmission_modes(2),'Location','northwest');
% set a title and key (legend) for the figure
title('Graph of Power Loss (6 thermal gen) against Power Demand');
%label axes
xlabel('Power Demand (MW)'); % x-axis label
ylabel('Power Loss (MW)');  % y-axis label
grid on;                    %enable grid on the graph
saveas(fig4,plot_Ploss);    % saves the graph to file

fig5 = figure('Visible', 'off');
y5 = [sw_loss;cond_loss;tl_loss]';
bar(demand,y5)
title('Switching, Conduction and Transmission Lossses for ED with HVDC'); %title
legend('Switching Losses','Conduction Losses','Transmission Losses','Location','northwest');
%label axes
xlabel('Loss Categories'); % x-axis label
ylabel('Power Loss (MW)'); % y-axis label
grid on;                    %enable grid on the graph
saveas(fig5,hvdc_Losses);    % saves the graph to file

fig6 = figure('Visible', 'off');
y6 = [F_cost(1,load_demand_values+1:load_demand_values*2);F_cost_inst(1,1:load_demand_values);F_cost(1,1:load_demand_values)]';
bar(demand,y6);
% set a title and key (legend) for the figure
title('Fuel Cost(6 thermal generators) against Power Demand'); %title
legend(strcat(transmission_modes(2),' stable'),strcat(transmission_modes(2),' with instability'),transmission_modes(1),'Location','northwest');
%label axes
xlabel('Power Demand (MW)'); % x-axis label
ylabel('Fuel Cost ($)'); % y-axis label
grid on;                    %enable grid on the graph
%legend(p(:,:),{'HVDC','HVAC'},'Location','north');
saveas(fig6,bar_Fcost);    % saves the graph to file

disp(['Success!!! Computed dispatch for ',num2str(load_demand_values)...
    , ' load demands. Results in this project folder with names: ',...
    num2str(doc_name),' , ',num2str(plot_Fcost),' , ',num2str(plot_Ploss)...
    ,' , ',num2str(plot_Iterr),' and ',num2str(hvdc_Losses)])

function[fuel_cost, generator_powers, P_loss]=eldnba(bestX)
global fuel_coefficients B power_demand Pg_limits transmission_type ...
    Cond_loss SW_loss TL_loss n f_cost beta tao time inst_const

%% Step 4a: Finding the power output(Pgi) from each generator
best_fit = abs(bestX);%Gives absolute value of element optimization_result

for i=1:n-1
    if best_fit(i)>1;
        best_fit(i)=1;
    else
    end
    P(i) =  Pg_limits(i+1,1) + best_fit(i)* (Pg_limits(i+1,2)-Pg_limits(i+1,1));
end
B11=B(1,1);
B1n=B(1,2:n);
Bnn=B(2:n,2:n);
A=B11;
BB1=2*B1n*P';
B1=BB1-1;
C1=P*Bnn*P';
C = power_demand - sum(P) + C1;
% finding and assigning roots of a polynomial to the varriable xl
x1= roots([A B1 C]);
% getting the absolute value of the minimum of the polynomial roots
Pgi = abs(min(x1));
% set Pgi to Pgi_max If value assigned to it exceeds maximum limit
if Pgi> sum(Pg_limits(:,2))
    Pgi = sum(Pg_limits(:,2));
elseif Pgi<Pg_limits(1,1)
    % set Pgi to Pgi_min If value assigned to it is below minimum limit
    Pgi=Pg_limits(1,1);
else
end
generator_powers=[Pgi P];

%% Step 4b: Finding the total power losses
TL_loss = 0.45*(generator_powers*B*generator_powers'); % transmission line loss
if strcmp(transmission_type,'HVDC') %assigning a different B for HVDC
    SW_loss = 0.2*(generator_powers*B*generator_powers'); % switching losses
    Cond_loss = 0.1*(generator_powers*B*generator_powers'); % conduction losses
    Conv_loss = SW_loss + Cond_loss; % conversion losses
    P_loss = TL_loss + Conv_loss; % total power loss HVDC
else
    P_loss = TL_loss; % total power loss HVAC
end
lam = abs(sum(generator_powers) - power_demand - P_loss);

%% Step 4c: Finding the total fuel cost
[F1,Ej,inst_const]= deal(zeros);
for i=1:n
    F1(i) = fuel_coefficients(i,1)*generator_powers(i)^2 +...
        fuel_coefficients(i,2)*generator_powers(i) + ...
        fuel_coefficients(i,3) + abs(fuel_coefficients(i,4)* ...
        sin(fuel_coefficients(i,5))*(Pg_limits(i,1) - generator_powers(i)));
end
if strcmp(transmission_type,'HVDC')
    for i=1:n
        Ej(i) = generator_powers(i)*time/60;
    end
    inst_const=sum(tao.*Ej)+sum(beta.*Ej);
end

fuel_cost = sum(F1) + 1000*lam;
f_cost = F1 + 1000*lam/n;
end

%% NOVEL BAT ALGORITHM
function [ bestX, bestIter ] = NBA( FitFunc, M, pop, dim, G, gamma, alpha, ...
    r0Max, r0Min, AMax, AMin, freqDMax, freqDMin, probMax, probMin, ...
    CMax, CMin, thetaMax, thetaMin, wMax, wMin )

% set the parameters
lb= -100 * ones( 1,dim );   % Lower bounds
ub= 100 * ones( 1,dim );    % Upper bounds
vLb = 0.6 * lb;
vUb = 0.6 * ub;

r = rand( pop, 1 ) .* 0.2 + 0;
r0 = rand( pop, 1 ) .* ( r0Max - r0Min ) + r0Min;
A = rand( pop, 1 ) .* ( AMax - AMin ) + AMin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialization

for i = 1 : pop
    x( i, : ) = lb + (ub - lb) .* rand( 1, dim );
    v( i, : ) = rand( 1, dim );
    fit( i ) = FitFunc( x( i, : ) );
end
pFit = fit; % The individual's best fitness value
pX = x;     % The individual's best position corresponding to the pFit

[ fMin, bestIndex ] = min( fit );  % fMin denotes the global optimum
% bestX denotes the position corresponding to fMin
bestX = x( bestIndex, : );
bestIter = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start the iteration.

for iteration = 1 : M
    
    % The compensation rates for Doppler effect in echoes
    C = rand( pop, 1 ) .* ( CMax - CMin ) + CMin;
    % The probability of habitat selection
    prob = rand( pop, 1 ) .* ( probMax - probMin ) + probMin;
    % ContractionCexpansion coefficient
    theta=( thetaMax - thetaMin ) * ( M - iteration )/(1.0 * M) + thetaMin;
    
    freqD = rand( pop, dim ) .* ( freqDMax - freqDMin ) + freqDMin;
    w = (wMax - wMin) * ( M - iteration )/(1.0 * M) + wMin; %Inertia weight
    
    meanP = mean( pX );
    meanA = mean( A );
    
    for i = 1 : pop
        if rand < prob
            if rand < 0.5
                x( i, : ) = bestX + theta * abs( meanP - pX(i, :) ) *...
                    log( 1.0/rand );
            else
                x( i, : ) = bestX - theta * abs( meanP - pX(i, :) ) *...
                    log( 1.0/rand );
            end
        else
            freqD( i, :) = freqD(i, :) .* ( 340 + v( i, : ) )./( 340 + ...
                v( bestIndex, : ) + realmin );
            v( i, : ) = w .* v( i, : ) + ( bestX - pX(i, :) ) .* ...
                freqD(i,:) .* ( 1 + C(i) .* ( bestX - pX(i, :) ) ./...
                ( abs( bestX - pX(i, :) ) + realmin ) );
            
            v( i, : ) = Bounds( v( i, : ), vLb, vUb );
            x( i, : ) = x( i, : ) + v( i, : );
        end
        
        % Local search
        if rand > r( i )
            randnValueA = randn( 1,dim ).* ( abs( A(i) - meanA )+ realmin);
            x( i, : ) = bestX .* ( 1 + randnValueA );
        end
        
        x( i, : ) = Bounds( x( i, : ), lb, ub );
        fit( i ) = FitFunc( x( i, : ) );
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the individual's best fitness vlaue and the global best one
    
    for i = 1 : pop
        if fit( i ) < pFit( i )
            pFit( i ) = fit( i );
            pX( i, : ) = x( i, : );
        end
        
        if( pFit( i ) < fMin && rand < A(i) )
            fMin = pFit( i );
            bestX = pX( i, : );
            bestIndex = i;
            bestIter = iteration;
            A(i) = A(i) * alpha;
            r(i) = r0(i) * ( 1 - exp( -gamma * iteration ) );
        end
    end
    
    if( iteration - bestIter > G )
        r = rand( pop, 1 ) .* 0.05 + 0.85;
        A = rand( pop, 1 ) .* ( AMax - AMin ) + AMin;
    end
    
end
% End of the main program

% The following functions are associated with the main program
% This function is the objective function
    function y = Sphere( x )
        y = sum( x .^ 2 );
    end
% Application of simple limits/bounds
    function s = Bounds( s, Lb, Ub)
        % Apply the lower bound vector
        temp = s;
        I = temp < Lb;
        temp(I) = Lb(I);
        
        % Apply the upper bound vector
        J = temp > Ub;
        temp(J) = Ub(J);
        % Update this new move
        s = temp;
    end
end