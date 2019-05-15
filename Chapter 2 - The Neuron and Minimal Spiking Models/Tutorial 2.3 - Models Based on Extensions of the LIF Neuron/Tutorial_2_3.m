%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Q1: LIF model with an adaptation current %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E_l = -75.0e-3;
V_th = -50.0e-3;
V_reset = -80.0e-3;

R_m = 100.0e6;

C_m = 100.0e-12;

E_k = -80.0e-3;

delta_G = 1.0e-9;

tau_sra = 0.2;

delta_t = 0.0001;

t = 0:delta_t:1.5;

V = zeros(1, length(t));
V(1) = E_l;

G_sra = zeros(1, length(t));

I_app = zeros(1, length(t));

I_app(5001:10001) = 500e-12;

%Use Forward Euler's method to integrate the equation
for i = 1:length(t)-1   
    if ( V(i) > V_th )         
        V(i) = V_reset;         
        G_sra(i) = G_sra(i) + delta_G;        
    end
    
    V(i+1) = V(i) + delta_t* ( (E_l-V(i))/R_m + G_sra(i)*(E_k - V(i)) + I_app(i))/C_m;

    G_sra(i+1) = G_sra(i) - delta_t*(G_sra(i)/tau_sra);
end

%Plot figure
figure(1)
subplot(3,1,1)
plot(t, 1e12*I_app);
ylabel('I_{app} (pA)');
subplot(3,1,2)
plot(t, 1000*V);
ylabel('V_m (mV)');
subplot(3,1,3)
plot(t, G_sra*1e9);
xlabel('Time (s)');
ylabel('G_{sra} (nS)');
figure();

% Part B: simulate model for 5s with a range of 20 different levels of 
% constant applied current (i.e. no step pulse)

t = 0:delta_t:5;

I_app = 240:5:550;

I_app = I_app * 1e-12;

%neuron_fires = zeros(1, length(I_app));

%firing_rate = zeros(1, length(I_app));

initialrate = zeros(size(I_app)); % array to store 1/(first ISI)
finalrate = zeros(size(I_app));   % array to store 1/(final ISI)
singlespike = zeros(size(I_app)); % array to store "1" for only 1 spike

%Loop through applied current values
for j = 1:length(I_app)   
    %Create vector to store membrane potential values
    V = zeros(1, length(t));
    V(1) = E_l;
    
    G_sra = zeros(1, length(t));
    
    spikes = zeros(size(t));
    
    %Loop through time vector
    for i = 1:length(t)-1
        if ( V(i) > V_th )         
            V(i) = V_reset;         
            G_sra(i) = G_sra(i) + delta_G;   
            %neuron_fires(j) = neuron_fires(j) + 1;
            spikes(i) = 1;
        end
    
        V(i+1) = V(i) + delta_t * ( (E_l-V(i))/R_m + G_sra(i)*(E_k - V(i)) + I_app(j))/C_m;

        G_sra(i+1) = G_sra(i) - delta_t*(G_sra(i)/tau_sra);
    end
    %Caclulate firing rate
    %firing_rate(j) = neuron_fires(j)/t(length(t));
    
     spiketimes = delta_t*find(spikes);           % extract the spike times
    
    if ( length(spiketimes) > 1 )           % if there is more than 1 spike
        ISIs = diff(spiketimes);            % ISI = interval between spikes
        initialrate(j) = 1/ISIs(1);     % inverse of first ISI
        if ( length(ISIs) > 1 )             % if there are further ISIs
            finalrate(j) = 1/ISIs(end); % inverse of final ISI
        end
        
    else
        if ( length(spiketimes) == 1 )      % if there is only one spike
            singlespike(j) = 1;         % record "1" for this trial
        end
    end
end

hold on;
plot(I_app*1e12, finalrate, 'k');

ISIindices = find(initialrate);
plot(1e12*I_app(ISIindices),initialrate(ISIindices),'ok');

ISIindices = find(singlespike);
plot(1e12*I_app(ISIindices),0*singlespike(ISIindices),'*k');

xlabel('I_{app} (nA)');
ylabel('Spike rate (Hz)');

legend('Final Rate',  '1/ISI(1)' , 'Single spike');

hold off;
figure();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Q2: AELIF model (Adaptive Exponential Leaky Integrate-and-Fire) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


