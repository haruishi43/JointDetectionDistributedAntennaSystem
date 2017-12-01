clear;
ccc_table = load('CCCtable_2antenna', ...
                 'CCCtable_conv_SINRp_alphap_QAMq_QAMp', ...   % no joint ml detection
                 'CCCtable_prop_SINRp_alphap_QAMq_QAMp');      % joint ml detection

%% Randomize:
rng('Shuffle');

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
noise_power = 1;
eirp = 0 + 30 - ( rnd + 10*log10( band ) );

% Scheduling parameters
num_select = 2;                     % # of user selected for each combination
[combination_table, tot_combinations] = create_combination_table( num_users, num_select );

%% Simulation parameters:
num_drops = 1;
time_interval = 10;
trial_per_drop = 1;

%% Initializing variables:
plr_from_bs_all = zeros(num_drops,num_cell,num_users);              % packet loss ratio 
channel_response_freq = zeros(num_users, num_cell, num_sc);         % channel response frequency
channel_response = zeros(num_users, num_cell, num_rb);

%% Create coordinates for each BS:
preset_coordinates = [1 3 6]; % For Coordinate Testing (has to chnage when num_users change)
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
        
        %% Round-Robin scheduling with Max-C   
        current_comb = 1;
        connection = 8 * ones(num_rb, num_users);
        
        signal = zeros(num_rb, num_select, num_select); % 1 is main signal, 2 is interference
        power = zeros(num_rb, num_select);
        alpha = zeros(num_rb, num_select);
        power_floor = zeros(num_rb, num_select);
        alpha_floor = zeros(num_rb, num_select);
        modulation = zeros(num_rb, num_select);
        modulation_jd = zeros(num_rb, num_select);
        ccc_output = zeros(num_rb, num_select);
        ccc_output_jd = zeros(num_rb, num_select);
        
        for rb = 1:num_rb
            % rr with max-c
            [ signal(rb, :, :), power(rb, :), connection(rb, :) ] = rr_max_c( num_users, combination_table(current_comb,:), all_signal_power(:, :, rb) );
               
            % calculate alpha and floor it (0 to 1 incremented by 0.1)
            [ alpha_floor(rb, :), alpha(rb, :) ] = calculate_alpha( squeeze( signal(rb, :, :) ) );
            
            % floor P/N (-10 to 30 incremented by 1)
            for i = 1:num_select
                power_floor(rb, i) = floor(power(rb, i));
                if power_floor(rb, i) >= 30
                    power_floor(rb, i) = 30;
                elseif power_floor(rb, i) <= -10
                    power_floor(rb, i) = -10;
                end
            end
            
            % find the best modulation
            modulation(rb, :) = find_best_mod( power_floor(rb, :), alpha_floor(rb, :), 0, ccc_table );
            modulation_jd(rb, :) = find_best_mod( power_floor(rb, :), alpha_floor(rb, :), 1, ccc_table );
            
            % look up CCC output
            ccc_output(rb, 1) = ccc_table.CCCtable_conv_SINRp_alphap_QAMq_QAMp( power_floor(rb, 1) + 11, round(10*(1-alpha_floor(rb, 1))) + 1, modulation(rb, 2), modulation(rb, 1));
            ccc_output(rb, 2) = ccc_table.CCCtable_conv_SINRp_alphap_QAMq_QAMp( power_floor(rb, 2) + 11, round(10*(1-alpha_floor(rb, 2))) + 1, modulation(rb, 1), modulation(rb, 2));
            
            ccc_output_jd(rb, 1) = ccc_table.CCCtable_prop_SINRp_alphap_QAMq_QAMp( power_floor(rb, 1) + 11, round(10*(1-alpha_floor(rb, 1))) + 1, modulation_jd(rb, 2), modulation_jd(rb, 1));
            ccc_output_jd(rb, 2) = ccc_table.CCCtable_prop_SINRp_alphap_QAMq_QAMp( power_floor(rb, 2) + 11, round(10*(1-alpha_floor(rb, 2))) + 1, modulation_jd(rb, 1), modulation_jd(rb, 2));
            
            % increment
            current_comb = current_comb + 1;
            if current_comb > tot_combinations
                current_comb = 1;
            end
        end
        toc
        
    end
    
end
