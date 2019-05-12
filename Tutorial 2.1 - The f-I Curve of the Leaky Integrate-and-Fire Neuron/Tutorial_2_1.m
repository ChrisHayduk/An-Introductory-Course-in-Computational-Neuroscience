%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Let's calculate membrane potential using the leaky integrate-and-fire neuron model%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Given values for the leaky integrate-and-fire neuron model
E_l = -70.0;
R_m = 5.0;
C_m = 2.0;

V_th = -50.0;
V_reset = -65.0;

%Create time vector
delta_t = 0.1;
t = 0:delta_t:100;

%Create vector to store membrane potentials at time t
%Set initial value of the membrane potential vector to E_l
V = zeros(1, length(t));
V(1) = E_l;

%Minimum applied current needed for the neuoron to produce spikes
I_th = 1/R_m * (V_th - E_l);

%Create vector of applied I (current) values
I_app = (I_th+0.1)*ones(length(t),1);

%Adding noise scaler
sigma_I = 1.0;

%Creating vector of noise variables
noise_vec = randn(size(t)) * sigma_I * sqrt(delta_t);

%Use Forward Euler's method to integrate the equation
for i = 2:length(t)
    dxdt = (E_l - V(i-1))/R_m + I_app(i);
    dxdt = dxdt*(1/C_m);
    
    V(i) = V(i-1) + delta_t * dxdt + noise_vec(i);
    
    if V(i) > V_th
        V(i) = V_reset;
    end
end

%Plot time vs. membrane potential
plot(t, V);
xlabel('Time');
ylabel('Membrane Potential');
figure()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Now let's calculate approximate firing rate vs. exact firing rate%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define new delta t
%Need smaller time steps in order to achieve firing rate of 100
delta_t = 0.001;
t = 0:delta_t:2;

%Create vector of I_app (applied current) values to be used
I_app = I_th:350:I_th+3200;

%Create vector to store number of neuron fires for each trial
neuron_fires = zeros(1, length(I_app));

%Create vector to store firing rate for each trial
firing_rate = zeros(1, length(I_app));

%Loop through applied current values
for j = 1:length(I_app)
    %Creating vector of noise variables
    noise_vec = randn(size(t)) * sigma_I * sqrt(delta_t);
    
    %Create vector to store membrane potential values
    V = zeros(1, length(t));
    V(1) = V_reset;
    %Loop through time vector
    for i = 2:length(t)
        dxdt = (E_l - V(i-1))/R_m + I_app(j);
        dxdt = dxdt*(1/C_m);
    
        V(i) = V(i-1) + delta_t * dxdt + noise_vec(i);
        
        %Check if membrane potential is above threshold
        if V(i) > V_th
           V(i) = V_reset;
           neuron_fires(j) = neuron_fires(j) + 1;
        end
    end
    
    %Caclulate firing rate
    firing_rate(j) = neuron_fires(j)/2;
end

%Create vector to store exact firing rate
firing_rate_exact = zeros(1, length(I_app));

%Create tau variable
tau_m = C_m * R_m;

%Loop through entries in firing rate exact vector
for i = 1:length(firing_rate_exact)
    %Calculate the two terms
    term1 = I_app(i) * R_m + E_l - V_reset;
    term2 = I_app(i) * R_m + E_l - V_th;
    
    %Check that the terms are greater than 0
    %So that log is defined
    if (term1 > 0) && (term2 > 0)
        %Generate final term
        final_term = (tau_m * log(term1) - tau_m * log(term2));
        
        %Check if final term is not infinity
        if final_term ~= Inf
            firing_rate_exact(i) = 1/final_term;
        end
    end
end


%Plot applied current vs. approximate firing rate
plot(I_app, firing_rate);
hold on;

%Plot applied current vs. exact firing rate
plot(I_app, firing_rate_exact);
xlabel('Applied Current');
ylabel('Firing Rate');
legend('Approximate','Exact')
hold off;
figure();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Now let's test how sigma effects our approximate firing rate%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create vector of sigma values
sigma_I = 0:10:100;

%Loop through sigma_I values
for k = 1:length(sigma_I)
    %Create vector to store number of neuron fires for each trial
    neuron_fires = zeros(1, length(I_app));

    %Create vector to store firing rate for each trial
    firing_rate = zeros(1, length(I_app));
    
    %Loop through applied current values
    for j = 1:length(I_app)
        %Creating vector of noise variables
        noise_vec = randn(size(t)) * sigma_I(k) * sqrt(delta_t);
    
        %Create vector to store membrane potential values
        V = zeros(1, length(t));
        V(1) = V_reset;
        %Loop through time vector
        for i = 2:length(t)
            dxdt = (E_l - V(i-1))/R_m + I_app(j);
            dxdt = dxdt*(1/C_m);
    
            V(i) = V(i-1) + delta_t * dxdt + noise_vec(i);
        
            %Check if membrane potential is above threshold
            if V(i) > V_th
                V(i) = V_reset;
                neuron_fires(j) = neuron_fires(j) + 1;
            end
        end
        
        %Caclulate firing rate
        firing_rate(j) = neuron_fires(j)/2;
    end 
    
    %Plot applied current vs. approximate firing rate with sigma included
    plot(I_app, firing_rate);
    hold on;
end

%Label
xlabel('Applied Current');
ylabel('Firing Rate');
legend(cellstr(num2str(sigma_I', 'Sigma = %-d')));
hold off;