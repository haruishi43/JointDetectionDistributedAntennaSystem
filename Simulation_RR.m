clear;
load('CCCtable_2antenna');

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
eirp = 0 + 30 - (rnd + 10*log10(band));

% Scheduling parameters
num_select = 2;                     % # of user selected for each combination
[combination_table, tot_combinations] = create_combination_table(num_users, num_select);

%% Simulation parameters:
num_drops = 1;
time_interval = 10;
trial_per_drop = 1;

%% Initializing variables:
plr_from_bs_all = zeros(num_drops,num_cell,num_users);              % packet loss ratio 
channel_response_freq = zeros(num_users, num_cell, num_sc);         % channel response frequency
channel_response = zeros(num_users, num_cell, num_rb);

%% Create coordinates for each BS:
preset_coordinates = [2 4 6]; % For Coordinate Testing
antenna_coordinates = create_bs_coordinate();

%% Simulation loop (change user placement):   
for drop = 1:num_drops
    
    %% Create Coordinates for each user:
    user_coordinates = create_user_coordinates( antenna_coordinates, 3, 500, preset_coordinates );
    
    %% Calculate Propagation Loss 
    plr_from_bs_all(drop, :, :) = create_plr_from_bs( antenna_coordinates, user_coordinates );
    
    %% Simulation loop (trial):
    for trial = 1:trial_per_drop
        
        %% Calculate Rayleigh Fading:
        channel_response_freq = add_rayleigh_fading( num_users, num_cell );
        
        %% Average to create channel response for each RB:
        all_signal_power = zeros(num_users, num_cell, num_rb);
        for user = 1:num_users
            for cell = 1:num_cell
                
                const = 10.^(( eirp  - plr_from_bs_all(drop, cell, user) ) / 10);
                
                for rb = 1:num_rb
                    
                    channel_response(user, cell, rb) = abs(mean(channel_response_freq(user, cell, num_sc_in_rb * (rb-1) + 1:num_sc_in_rb * rb)));

                    % signal in real number domain
                    all_signal_power(user, cell, rb) = sqrt(shadowing_var)*10^(randn(1,1)) * const * ( abs(channel_response(user, cell, rb)).^2 );

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
        ccc_output = zeros(num_rb, num_select);
        
        for rb = 1:num_rb
            % signal of the rb
            cs = all_signal_power(:, :, rb);
            
            cc = combination_table(current_comb,:);
            
            % first user
            user1 = cc(1);
            user2 = cc(2);
            
            % max-c (best signal power is chosen)
            cs_user = cs(user1, : );
            [s, index] = max( cs_user );
            
            connection(rb, user1) = index;
            signal(rb, 1, 1) = s;
            signal(rb, 2, 2) = cs(user2, index);
            
            % add to power (to dB)
            % noise is already calculated through rnd
            power(rb, 1) = 10*log10(abs(s));
            power(rb, 2) = 10*log10(abs(signal(rb, 2, 2)));
            
            cs_user = cs(user2, :);
            [s, index] = max( cs_user );
            
            % find if the DA is already being used or not
            used = find(index == connection(rb, :));
            if isempty(used) == 1
                connection(rb, user2) = index;
                
                signal(rb, 2, 1) = s;
                signal(rb, 1, 2) = cs(user1, index);
                
                power(rb, 2) = power(rb, 2) + 10*log10(abs(s));
                power(rb, 1) = power(rb, 1) + 10*log10(abs(signal(rb, 1, 2)));
            end
               
            
            % calculate alpha
            alpha1 = signal(rb, 1, 1) / (signal(rb, 1, 1) + signal(rb, 1, 2));    
            alpha2 = signal(rb, 2, 1) / (signal(rb, 2, 1) + signal(rb, 2, 2));
            
            if isnan(alpha1)
                alpha1 = 1;
            end
            if isnan(alpha2)
                alpha2 = 1;
            end
            
            alpha(rb, :) = [alpha1 alpha2];
            
            
            % floor P/N and alpha
            for i = 1:num_select
                power_floor(rb, i) = floor(power(rb, i));
                if power_floor(rb, i) >= 30
                    power_floor(rb, i) = 30;
                elseif power_floor(rb, i) <= -10
                    power_floor(rb, i) = -10;
                end
            end
            
            % round down
            alpha_floor(rb, :) = floor(10*alpha(rb, :)) / 10;
            
            % find the best modulation
            mod_list1 = squeeze(CCCtable_conv_SINRp_alphap_QAMq_QAMp( power_floor(rb, 1) + 11, 10*alpha_floor(rb, 1) + 1, :, :));
            mod_list2 = squeeze(CCCtable_conv_SINRp_alphap_QAMq_QAMp( power_floor(rb, 2) + 11, 10*alpha_floor(rb, 2) + 1, :, :));
            
            tot_mod_list = mod_list1 + mod_list2';
            [value, index] = max(tot_mod_list(:));
            [row, col] = ind2sub(size(tot_mod_list), index);
            
            modulation(rb, :) = [col, row];
            
            ccc_output(rb, 1) = CCCtable_conv_SINRp_alphap_QAMq_QAMp( power_floor(rb) + 11, 10*alpha_floor(rb, 1) + 1, modulation(rb, 2), modulation(rb, 1));
            ccc_output(rb, 2) = CCCtable_conv_SINRp_alphap_QAMq_QAMp( power_floor(rb) + 11, 10*alpha_floor(rb, 2) + 1, modulation(rb, 1), modulation(rb, 2));
            
            % increment
            current_comb = current_comb + 1;
            if current_comb > tot_combinations
                current_comb = 1;
            end
        end
        
        
        
        
    end
    
end
