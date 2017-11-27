clear;

%% Randomize:
%rng('Shuffle');

%% Model parameters:
num_users = 3;                      % # of users
num_cell = 7;                       % # of cell

num_rb = 24;                        % # of resource blocks in 1 OFDM symbol
num_sc_in_rb = 12;                  % # of subcarriers in resource blocks
num_sc = num_rb * num_sc_in_rb;     % # of total subcarriers

band_per_rb = 180*10^3;             % frequency band range for each rb (Hz)
band = band_per_rb * num_rb;        % total frequency band

shadowing_ave = 0;
shadowing_var = 8;
rnd = -174;                         % Reciever Noise Density

%% Simulation parameters:
num_drops = 1;
time_interval = 10;
trial_per_drop = 1;
saving = 1;

%% Initializing variables:
plr_from_bs_all = zeros(num_drops,num_cell,num_users);              % packet loss ratio 
channel_response_freq = zeros(num_users, num_cell, num_sc);         % channel response frequency
channel_response = zeros(num_users, num_cell, num_rb);

%% Create coordinates for each BS:
preset_coordinates = [1 2 7]; % For Coordinate Testing
antenna_coordinates = create_bs_coordinate();

%% Simulation loop (change user placement):   
for drop = 1:num_drops
    
    %% Create Coordinates for each user:
    user_coordinates = create_user_coordinates( antenna_coordinates, 3, 500, preset_coordinates );
    
    %% Calculate Packet Loss Ratio
    plr_from_bs_all(drop, :, :) = create_plr_from_bs( antenna_coordinates, user_coordinates );
    
    %% Simulation loop (trial):
    for trial = 1:trial_per_drop
        
        %% Calculate Rayleigh Fading:
        channel_response_freq = add_rayleigh_fading( num_users, num_cell );
        
        %% Average to create channel response for each RB:
        for user = 1:num_users
            for cell = 1:num_cell
               for rb = 1:num_rb
                    channel_response(user, cell, rb) = abs(mean(channel_response_freq(user, cell, num_sc_in_rb * (rb-1) + 1:num_sc_in_rb * rb)));
               end
            end
        end 
        
        
        
        
        
    end
    
end
