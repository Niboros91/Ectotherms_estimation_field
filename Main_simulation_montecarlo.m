% ODE model for Drosophila suzukii
%Initialization
clear, clc % To clear the workspace and the command function

Montecarlo_iterations = 100; %100

%% Traps characteristics

% In between 5 and 10 days
 trap_time = 7;
 trap_on = 0; %Initially set to 0
 trap_attrac = 0.15; %Ratio of attracted insects from the trap
%% PARAMETERS functions

%Growth rate (Drosophila suzukii)
 a = 1.2*(10^(-4));
 T_l = 3;
 T_m = 30;
 m = 6;
 
%Uncertainty growth rate
v_a = 0.15 * (10^(-4));
v_T_l = 2;
v_T_m = 1;
v_m = 3;

%Mortality rate (Drosophila suzukii)

a1 = -5.4E-06;
b1 = 0.0005194;
c1 = -0.0116827;
d1 = 2.16E-05;
e1 = 1.3146586;

v_a1 = 1.4E-05;
v_b1 = 0.0008;
v_c1 = 0.02;
v_d1 = 0.2;
v_e1 = 0.9;

% Birth rate (Drosophila suzukii)
alpha = 659.06;
gamma = 88.53;
lambda = 52.32;
delta = 6.06;
tau = 22.87;

% Sex ratio
S_R = 0.5;

% Mating ratio
R_remate = 0;
R_mate = 1;

% Additional parameters
N_stages = 9; % Eggs/3 larva stages/ pupa /male/ unmated female/ mated female/ traps

    
    %% DATA
  
    %Load additional data. To update depending on the temperature used

    Field = 'Montelibretti'; %Montelibretti Monterotondo
    Year = '2018'; % 2017 2018 2019
    Folder = 'data_field';
    time_diff = 0;
    
    if contains(Year,'2019')
        if contains (Field,'Montelibretti')
            time_diff = 35;
        else
            time_diff = 37;
        end
    end
    
    temperature = strcat(pwd,'\',Folder,'\','Temperatures_',Field,'_',Year,'.csv');
    population_data = strcat(pwd,'\',Folder,'\','Adults_',Field,'_',Year,'.csv');
    Temp_avg = csvread(temperature, 1);
    valid_1 = csvread(population_data, 1);
    
    
    %Separe measurements
    measures= valid_1;
    traps = measures(:,1);
    counter_trap =1;

    %% MEASUREMENT MODEL

    Trap_eff = 1;

    % We assume that we can measure only the adult male states
    
   C = [zeros(1,5), 0, 0, 0, Trap_eff];
    
   %% ERRORS STORAGE
   
   error_total = zeros(length(measures)-1, length(measures)-1);
   error_prediction = zeros(1, length(measures)-1);
   error_open = zeros(1, length(measures)-1);
   

    
    %% SIMULATIONS

variation = 0.20; %Variation of 10% up and down

error_prediction_hist = zeros(Montecarlo_iterations,(length(measures)-1)); 
error_open_hist = zeros(Montecarlo_iterations,(length(measures)-1));

counter_test = 0;
x_open_hist_all = [];
x_ekf_hist_all = [];

for it_mont = 1:Montecarlo_iterations    
   % Adding the variation to each iteration
    x_ini = initial_conditions(Field, Year,N_stages);
    x_ini_open = initial_conditions(Field, Year,8);

    egg_var = (-x_ini(1)*variation) + ((x_ini(1)*variation)-(-x_ini(1)*variation)) .* rand(1,1);
    adult_var = (-x_ini(8)*variation) + ((x_ini(8)*variation)-(-x_ini(8)*variation)) .* rand(1,1);
    
    x_ini(1) = x_ini(1) + egg_var;
    x_ini(8) = x_ini(8) + adult_var;
    x_ini_open(1) = x_ini_open(1) + egg_var;
    x_ini_open(8) = x_ini_open(8) + adult_var;
    
    

        %Initialize the counters
        counter =1; %counter of the measurements available  
        counter_meas = 1; %counter of the measurement compared
        counter_check = 1;


        % REAL RATES
        %Compute the rates based on temperature
        G_R = growth_rate(Temp_avg(1),a,T_l,T_m,m); %Calling the growth rate function

        B_R = birth_rate_suzuki(Temp_avg(1),alpha,gamma,lambda,delta,tau);%Calling the birth rate function

        M_R = mortality_rate(Temp_avg(1),a1,b1,c1,d1,e1); %Calling the mortality rate function

        %Initialize stages (synthetic data)
        Pest_stages(1:N_stages) = stage; % We create an array of the stage class
        Pest_stages = Initialize_stages_ode(B_R,M_R,G_R,S_R,R_mate,R_remate,Pest_stages); %We initialize the parameters associated to each stage class

        % MEASURED uncertainty rates  
        [w_d, w_m, w_f] = rate_noise(Temp_avg(1),v_a1,v_b1,v_c1,v_d1,v_e1,v_a,v_T_l,v_T_m,v_m,a,T_l,T_m,m);

        uncertainty_rate = [w_d^2, w_m^2, w_f^2, w_d^2, w_m^2, w_d^2, w_m^2, w_d^2, w_m^2, w_d^2, w_m^2, w_d^2, w_m^2, w_d^2, w_m^2];

        %Update the matrix A (synthetic data)
        A_cont =compute_A_continous_integral(Pest_stages,trap_on);
        sysc = ss(A_cont,[],eye(9),[]);
        sysd = c2d(sysc,1,'zoh');
        A_dis = sysd.A;
        
        %Update the matrix A (synthetic data)
        A_cont_open =compute_A_continous(Pest_stages);
        sysc_open = ss(A_cont_open,[],eye(8),[]);
        sysd_open = c2d(sysc_open,1,'zoh');
        A_dis_open = sysd_open.A;


        %Initial estimation (Open loop and EKF)
        x =x_ini;
        x_open = x_ini_open;

        % Historic evolution of the states
        x_hist= x; %We start the array x_hist to keep a historic of the evolution of the states
        x_open_hist = x_open;

        %Historic evolution of the covariance
        P =  eye(N_stages)*1;

        % Jacobian state
   L = [-x(1) -x(1) x(8) zeros(1,13) ;... %Egg (1)
    x(1) 0 0 -x(2) -x(2) zeros(1,11);... %L1 (2)
    zeros(1,3) x(2) 0 -x(3) -x(3) zeros(1,9);... %L2 (3)
    zeros(1,5) x(3) 0 -x(4) -x(4) zeros(1,7);... %L3 (4)
    zeros(1,7) x(4) 0 -x(5) -x(5) zeros(1,5);...%P (5)
    zeros(1,9) 0.5*x(5) 0 -x(6) zeros(1,4);...%AM (6)
    zeros(1,9) 0.5*x(5) 0 0  -x(7) -x(7) 0 0;%NMF (7)
    zeros(1,12) x(7) 0 -x(8) 0; %MF (8)
    zeros(1,15) 1];

        Simulation_time =length(Temp_avg); %Simulation lenght based on the temperature array introduced
        for t=1:(Simulation_time-1) %Loop for the simulation

            x_open = A_dis_open*x_open; %Compute the state at the next time step

            measuring = 0;
            y=0;
            
            if measures(counter)-trap_time < t
                trap_on = trap_attrac;            
            end
            if counter <= length(measures)
                if measures(counter,1)-time_diff == t
                    y = measures(counter,2);
                    measuring = 1;
                    counter = counter +1;
                end
            end
        
            [x,P] = EKF_nonlinear_observation_integral(t,A_dis,L,C,x_ini,y,measuring,P);
            
            if measuring ==1 %We rest the traps
               x(9) = 0; 
               trap_on = 0;
            end

            x_hist = [x_hist,x]; %Store the states
            x_open_hist = [x_open_hist,x_open]; %Store the states

            % Recompute rates for the new conditions
            G_R = growth_rate(Temp_avg(t+1),a,T_l,T_m,m);
            B_R = birth_rate_suzuki(Temp_avg(1),alpha,gamma,lambda,delta,tau);
            M_R = mortality_rate(Temp_avg(t+1),a1,b1,c1,d1,e1);
            Pest_stages = Initialize_stages_ode(B_R,M_R,G_R,S_R,R_mate,R_remate,Pest_stages);


            % Rates uncertainty based on current conditions
             [w_d, w_m, w_f] = rate_noise(Temp_avg(t+1),v_a1,v_b1,v_c1,v_d1,v_e1,v_a,v_T_l,v_T_m,v_m,a,T_l,T_m,m);

             uncertainty_rate = [w_d, w_m, w_f, w_d, w_m, w_d, w_m, w_d, w_m, w_d, w_m, w_d, w_m, w_d, w_m,];


            %Update the matrix A (synthetic data)
            A_cont =compute_A_continous_integral(Pest_stages,trap_on);
            sysc = ss(A_cont,[],eye(9),[]);
            sysd = c2d(sysc,1,'zoh');
            A_dis = sysd.A;

            %Update the matrix A (open loop)
            A_cont_open =compute_A_continous(Pest_stages);
            sysc_open = ss(A_cont_open,[],eye(8),[]);
            sysd_open = c2d(sysc_open,1,'zoh');
            A_dis_open = sysd_open.A;
        
        L = [-x(1) -x(1) x(8) zeros(1,13) ;... %Egg (1)
            x(1) 0 0 -x(2) -x(2) zeros(1,11);... %L1 (2)
            zeros(1,3) x(2) 0 -x(3) -x(3) zeros(1,9);... %L2 (3)
            zeros(1,5) x(3) 0 -x(4) -x(4) zeros(1,7);... %L3 (4)
            zeros(1,7) x(4) 0 -x(5) -x(5) zeros(1,5);...%P (5)
            zeros(1,9) 0.5*x(5) 0 -x(6) zeros(1,4);...%AM (6)
            zeros(1,9) 0.5*x(5) 0 0  -x(7) -x(7) 0 0;%NMF (7)
            zeros(1,12) x(7) 0 -x(8) 0; %MF (8)
            zeros(1,15) 1];

        end
    
        x_open_hist_all = [x_open_hist_all;x_open_hist(6,:)];
        x_ekf_hist_all = [x_ekf_hist_all;x_hist(6,:)];
    end 


%% PLOTS
t1 = datetime(2020,3,15,12,0,0);
t=t1+days(0:Simulation_time-1);


tt = 1:length(Temp_avg);
t_measures = 1:(length(measures)-1);

title_plot = strcat(Field,32,Year);

% Plot adult males population  

figure
plot(t,x_hist(9,:),'LineWidth',2);
hold on
scatter(t(valid_1(:,1)-time_diff),valid_1(:,2),'LineWidth',2)
plot(t(valid_1(:,1)-time_diff),valid_1(:,2),":",'LineWidth',1.5)
legend("EKF traps","Field data",'Fontsize',15);
xlabel('Time [days]','Fontsize',15);
ylabel('Number of adult males','Fontsize',15);
title(title_plot,'Fontsize',20);
set(gca,'FontSize',15)

%%

figure
plot(t,x_open_hist_all(:,:),'LineWidth',2);
hold on
scatter(t(valid_1(:,1)-time_diff),valid_1(:,2),'LineWidth',2)
plot(t(valid_1(:,1)-time_diff),valid_1(:,2),":",'LineWidth',1.5)
xlabel('Time [days]','Fontsize',15);
ylabel('Number of adult males','Fontsize',15);
title('Open loop','Fontsize',20);
set(gca,'FontSize',15)

figure
plot(t,x_ekf_hist_all(:,:),'LineWidth',2);
hold on
scatter(t(valid_1(:,1)-time_diff),valid_1(:,2),'LineWidth',2)
plot(t(valid_1(:,1)-time_diff),valid_1(:,2),":",'LineWidth',1.5)
xlabel('Time [days]','Fontsize',15);
ylabel('Number of adult males','Fontsize',15);
title('EKF','Fontsize',20);
set(gca,'FontSize',15)
