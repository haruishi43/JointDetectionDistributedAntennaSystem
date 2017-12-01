clear;
ccc_table = load('CCCtable_2antenna', ...
                 'CCCtable_conv_SINRp_alphap_QAMq_QAMp', ...   % no joint ml detection
                 'CCCtable_prop_SINRp_alphap_QAMq_QAMp');      % joint ml detection

%% Randomize:
rng('Shuffle');

%% Model parameters:
num_users = 3;                      % # of users
num_cell = 7;                       % # of cell
preset_coordinates = [2 4 6];   % For Coordinate Testing (has to change when num_users change)

num_rb = 24;                        % # of resource blocks in 1 OFDM symbol
num_sc_in_rb = 12;                  % # of subcarriers in resource blocks
num_sc = num_rb * num_sc_in_rb;     % # of total subcarriers

band_per_rb = 180*10^3;             % frequency band range for each rb (Hz)
band = band_per_rb * num_rb;        % total frequency band

shadowing_ave = 0;
shadowing_var = 8;
rnd = -174;                         % Reciever Noise Density
noise_power = 1;
eirp = 0 + 30 - ( rnd + 10*log10( band ) );

% Scheduling parameters
num_select = 2;                     % # of user selected for each combination
[combination_table, tot_combinations] = create_combination_table( num_users, num_select );

%% Simulation parameters:
num_drops = 1;
time_interval = 1;
trial_per_drop = 1;

%% Initializing variables:
plr_from_bs_all = zeros(num_drops,num_cell,num_users);              % packet loss ratio 
channel_response_freq = zeros(num_users, num_cell, num_sc);         % channel response frequency
channel_response = zeros(num_users, num_cell, num_rb);

%% Create coordinates for each BS:
antenna_coordinates = create_bs_coordinate();

%% Simulation loop (change user placement):   
for drop = 1:num_drops
    
    %% Create Coordinates for each user:
    user_coordinates = create_user_coordinates( antenna_coordinates, num_users, 500, preset_coordinates );
    
    %% Calculate Propagation Loss 
    plr_from_bs_all(drop, :, :) = create_plr_from_bs( antenna_coordinates, user_coordinates );
    
    %% Simulation loop (trial):
    for trial = 1:trial_per_drop
        tic
        
        %% Calculate Rayleigh Fading:
        channel_response_freq = add_rayleigh_fading( num_users, num_cell );
        
        %% Average to create channel response for each RB:
        all_signal_power = zeros(num_users, num_cell, num_rb);
        for user = 1:num_users
            for cell = 1:num_cell
                
                const = 10.^(( eirp  - plr_from_bs_all(drop, cell, user) ) / 10);
                
                for rb = 1:num_rb
                    
                    % channel response (average of all subcarriers in a
                    % resource block
                    channel_response(user, cell, rb) = mean( channel_response_freq( user, cell, num_sc_in_rb * (rb-1) + 1:num_sc_in_rb * rb ) );

                    % signal in real number domain
                    all_signal_power(user, cell, rb) = sqrt(shadowing_var)*10^( randn(1,1) ) * const * ( abs( channel_response(user, cell, rb) ).^2 );

                end
            end
        end 
        
         
        current_comb = 1;   % for incrementing
        connection = 8 * ones(time_interval, num_rb, num_users);
        ccc_output = zeros(time_interval, num_rb, num_select);
        ccc_output_jd = zeros(time_interval, num_rb, num_select);
        
        for t = 1:time_interval
            power_floor = zeros(num_rb, num_select);
            alpha_floor = zeros(num_rb, num_select);
            
            for rb = 1:num_rb
              %% Round-Robin scheduling with Max-C  
                [ ccc_output(t, rb, :), ccc_output_jd(t, rb, :), power_floor(rb, :), alpha_floor(rb, :), connection(t, rb, :) ] = rr_max_c( num_users, combination_table(current_comb,:), all_signal_power(:, :, rb), ccc_table );
                
              %% Round-Robin schedulig with Max-C/I
                
                
                % increment
                current_comb = current_comb + 1;
                if current_comb > tot_combinations
                    current_comb = 1;
                end
            end
        end
        toc
        sum(sum(sum(ccc_output)))
        sum(sum(sum(ccc_output_jd)))
    end
    
end
