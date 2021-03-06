%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Begin by initializing all constants %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Given values for the leaky integrate-and-fire neuron model
E_l = -70.0;
R_m = 100.0;
C_m = 0.1;

V_th = -50.0;
V_reset = -65.0;

%Create time vector
delta_t = 0.00001;
t = 0:delta_t:2;

%Create vector of applied I (current) values
I_app = 100:100:600;

% Q1 - Forced Voltage Clamp

%Create vector to store number of neuron fires for each trial
neuron_fires = zeros(1, length(I_app));

%Create vector to store firing rate for each trial
firing_rate_Q1 = zeros(1, length(I_app));

%Noise variable
sigma_I = 0.0;

%Length of voltage clamp in milliseconds
tau_ref = 0.0025;

%Create vector to store mean membrane potential for each trial
mean_membrane_potential_Q1 = zeros(1, length(I_app));

%Loop through applied current values
for j = 1:length(I_app)
    %Variable to store time of last spike
    time_of_last_spike = -tau_ref - 1;

    %Creating vector of noise variables
    noise_vec = randn(size(t)) * sigma_I * sqrt(delta_t);
    
    %Create vector to store membrane potential values
    V = zeros(1, length(t));
    V(1) = E_l;
    
    %Loop through time vector
    for i = 2:length(t)
        dxdt = (E_l - V(i-1))/R_m + I_app(j);
        dxdt = dxdt*(1/C_m);
    
        V(i) = V(i-1) + delta_t * dxdt + noise_vec(i);
        
        %Check if neuron is still within the refractory period
        if t(i) < time_of_last_spike + tau_ref
            V(i) = V_reset;
        end
        
        %Check if membrane potential is above threshold
        if V(i) > V_th
           time_of_last_spike = t(i);
           V(i) = V_reset;
           neuron_fires(j) = neuron_fires(j) + 1;
        end
        
    end
    
    %Caclulate firing rate
    firing_rate_Q1(j) = neuron_fires(j)/2;
    
    mean_membrane_potential_Q1(j) = mean(V);
end


% Q2 - Threshold Increase


tau_vth = 0.001;

neuron_fires = zeros(1, length(I_app));

%Create vector to store firing rate for each trial
firing_rate_Q2 = zeros(1, length(I_app));

%Create vector to store mean membrane potential for each trial
mean_membrane_potential_Q2 = zeros(1, length(I_app));

%Loop through applied current values
for j = 1:length(I_app)
    %Creating vector of noise variables
    noise_vec = randn(size(t)) * sigma_I * sqrt(delta_t);
    
    %Create vector to store membrane potential values
    V = zeros(1, length(t));
    V(1) = E_l;
    
    V_th = ones(1, length(t))*-50;
    
    %Loop through time vector
    for i = 2:length(t)
        dxdt = (E_l - V(i-1))/R_m + I_app(j);
        dxdt = dxdt*(1/C_m);
    
        V(i) = V(i-1) + delta_t * dxdt + noise_vec(i);
        
        V_th(i) = V_th(i-1) + delta_t * (-50.0 - V_th(i-1))/tau_vth;
                
        %Check if membrane potential is above threshold
        if V(i) > V_th(i)
           V(i) = V_reset;
           V_th(i) = 200;
           neuron_fires(j) = neuron_fires(j) + 1;
        end
        
    end
    
    %Caclulate firing rate
    firing_rate_Q2(j) = neuron_fires(j)/2;
    
    mean_membrane_potential_Q2(j) = mean(V);
end


% Q3 - Refactory Conductance with Threshold Increase

E_k = -80.0;

tau_gref = 0.0002;

neuron_fires = zeros(1, length(I_app));

%Create vector to store firing rate for each trial
firing_rate_Q3 = zeros(1, length(I_app));

%Create vector to store mean membrane potential for each trial
mean_membrane_potential_Q3 = zeros(1, length(I_app));

%Loop through applied current values
for j = 1:length(I_app)
    %Creating vector of noise variables
    noise_vec = randn(size(t)) * sigma_I * sqrt(delta_t);
    
    %Vector to store refactory conductance
    G_ref = zeros(1, length(t));
    
    %Create vector to store membrane potential values
    V = zeros(1, length(t));
    V(1) = E_l;
    
    V_th = ones(1, length(t))*-50;
    
    %Loop through time vector
    for i = 2:length(t)
        dxdt = (E_l - V(i-1))/R_m + (G_ref(i-1) * (E_k - V(i-1))) + I_app(j);
        dxdt = dxdt*(1/C_m);
            
        V(i) = V(i-1) + delta_t * dxdt + noise_vec(i);
        
        G_ref(i) = G_ref(i-1) - (delta_t * G_ref(i-1)/tau_gref);
        
       
        %Check if membrane potential is above threshold
        if V(i) > V_th(i)
           V(i) = V_reset;
           V_th(i) = 200;
           G_ref(i) = G_ref(i) + 2;
           neuron_fires(j) = neuron_fires(j) + 1;
        end 
        
    end
    
    %Caclulate firing rate
    firing_rate_Q3(j) = neuron_fires(j)/2;
    
    mean_membrane_potential_Q3(j) = mean(V);
    
end

%Plot applied current vs. firing rate for all 3 models
plot(I_app, firing_rate_Q1);
hold on;
plot(I_app, firing_rate_Q2);
plot(I_app, firing_rate_Q3);
xlabel('Applied Current');
ylabel('Firing Rate');
legend('Voltage Clamp','Threshold Increase', 'Refactory Conductance and Threshold Increase')
hold off;
figure();

%Plot applied current vs. mean membrane potential for all 3 models
plot(I_app, mean_membrane_potential_Q1);
hold on;
plot(I_app, mean_membrane_potential_Q2);
plot(I_app, mean_membrane_potential_Q3);
xlabel('Applied Current');
ylabel('Mean Membrane Potential');
legend('Voltage Clamp','Threshold Increase', 'Refactory Conductance and Threshold Increase')
hold off;
figure();

%Plot firing rate vs. mean membrane potential for all 3 models
plot(firing_rate_Q1, mean_membrane_potential_Q1);
hold on;
plot(firing_rate_Q2, mean_membrane_potential_Q2);
plot(firing_rate_Q3, mean_membrane_potential_Q3);
xlabel('Firing Rate');
ylabel('Mean Membrane Potential');
legend('Voltage Clamp','Threshold Increase', 'Refactory Conductance and Threshold Increase')
hold off;