clear;
ccc_table = load('CCCtable_2antenna', ...
                 'CCCtable_conv_SINRp_alphap_QAMq_QAMp', ...   % no joint ml detection
                 'CCCtable_prop_SINRp_alphap_QAMq_QAMp');      % joint ml detection

%% Randomize:
rng('Shuffle');

%% Model parameters:
num_users = 10;                      % # of users
num_cell = 7;                       % # of cell
num_outer_macro = 6;
preset_coordinates = [1 3 5 6 7];   % For Coordinate Testing (has to change when num_users change)

num_rb = 24;                        % # of resource blocks in 1 OFDM symbol
num_sc_in_rb = 12;                  % # of subcarriers in resource blocks
num_sc = num_rb * num_sc_in_rb;     % # of total subcarriers

band_per_rb = 180*10^3;             % frequency band range for each rb (Hz)
band = band_per_rb * num_rb;        % total frequency band

shadowing_ave = 0;
shadowing_var = 8;
rnd = -174;                         % Reciever Noise Density
noise_power = rnd + 10*log10( band );
eirp = 0 + 30;

% Scheduling parameters
num_select = 2;                     % # of user selected for each combination
[combination_table, tot_combinations] = create_combination_table( num_users, num_select );

%% Simulation parameters:
num_drops = 5;
time_interval = 10;
trial_per_drop = 5;

%% Initializing variables:
plr_from_bs_all = zeros(num_drops, num_users, num_cell);              % propagation loss ratio
plr_from_outer_cell = zeros(num_drops, num_outer_macro, num_users, num_cell);          % propagation loss ratio from outer cell

channel_response_freq = zeros(num_users, num_cell, num_sc);           % channel response frequency
channel_response = zeros(num_users, num_cell, num_rb);

%% Saving variables:
throughput_one_user = zeros(num_drops, trial_per_drop, time_interval, num_rb);
% throughput = zeros(num_drops, trial_per_drop, time_interval, num_rb, num_select);
% throughput_jd = zeros(num_drops, trial_per_drop, time_interval, num_rb, num_select);
throughput_ci = zeros(num_drops, trial_per_drop, time_interval, num_rb, num_select);
throughput_ci_jd = zeros(num_drops, trial_per_drop, time_interval, num_rb, num_select);

%% Create coordinates for each BS:
antenna_coordinates = create_bs_coordinate();

%% Create outer cell coordinates:
outer_cell_coordinates = create_outer_cell_coordinates();

%% Simulation loop (change user placement):   
for drop = 1:num_drops
    
    %% Create Coordinates for each user:
    user_coordinates = create_user_coordinates( antenna_coordinates, num_users );
    
    %% Calculate Propagation Loss 
    plr_from_bs_all(drop, :, :) = create_plr_from_bs( antenna_coordinates, user_coordinates );
    plr_from_outer_cell(drop, :, :, :) = create_plr_from_outer_cell( outer_cell_coordinates, user_coordinates );
    
    %% Simulation loop (trial):
    for trial = 1:trial_per_drop
        tic
        
        %% Calculate Rayleigh Fading:
        channel_response_freq = add_rician_fading( num_users, num_cell );
        
        %% Average to create channel response for each RB:
        all_signal_power = zeros(num_users, num_cell, num_rb);
        for user = 1:num_users
            for cell = 1:num_cell
                
                const = 10.^(( eirp  - plr_from_bs_all(drop, user, cell) ) / 10);
                
                for rb = 1:num_rb
                    
                    % channel response (average of all subcarriers in a
                    % resource block) 
                    channel_response(user, cell, rb) = mean( channel_response_freq( user, cell, num_sc_in_rb * (rb-1) + 1:num_sc_in_rb * rb ) );
                    % this might be better:
                    %channel_response(user, cell, rb) = mean( abs( channel_response_freq( user, cell, num_sc_in_rb * (rb-1) + 1:num_sc_in_rb * rb ) ).^2 );
                    
                    % signal in real number domain
                    all_signal_power(user, cell, rb) = 10^( sqrt(shadowing_var)*randn(1,1) / 10 ) * const * ( abs( channel_response(user, cell, rb) ).^2 );
                    
                end
            end
        end
        
        %% Signal power from outer cell:
        all_signal_power_outer = zeros(num_outer_macro, num_users, num_cell, num_rb);
        for macro = 1:num_outer_macro
            channel_response_macro = zeros(num_users, num_cell, num_rb);
            channel_response_macro_freq = add_rayleigh_fading( num_users, num_cell );
            
            for user = 1:num_users
                for cell = 1:num_cell
                    
                    const = 10.^(( eirp  - plr_from_outer_cell(drop, macro, user, cell) ) / 10);
                    
                    for rb = 1:num_rb
                        
                        channel_response_macro(user, cell, rb) = mean( channel_response_macro_freq( user, cell, num_sc_in_rb * (rb-1) + 1:num_sc_in_rb * rb ) );
                        % signal in real number domain
                        all_signal_power_outer(macro, user, cell, rb) = const * ( abs( channel_response_macro(user, cell, rb) ).^2 );

                    end
                end
            end
        end
        
        
        %% Scheduling:
        current_user = 1;   % for incrementing single user (start from user 1)
        current_comb = 1;   % for incrementing combination 
        
        connection = 8 * ones(time_interval, num_rb, num_users);
        connection_jd = 8 * ones(time_interval, num_rb, num_users);
        
        ccc_output_one_user = zeros(time_interval, num_rb);
%         ccc_output = zeros(time_interval, num_rb, num_select);
%         ccc_output_jd = zeros(time_interval, num_rb, num_select);
        ccc_output_ci = zeros(time_interval, num_rb, num_select);
        ccc_output_ci_jd = zeros(time_interval, num_rb, num_select);
        
        for t = 1:time_interval
            for rb = 1:num_rb
              %% Max-C(/I) scheduling for single user
                ccc_output_one_user(t, rb) = single_user_scheduling( current_user, all_signal_power(:, :, rb), all_signal_power_outer(:, :, :, rb), noise_power, ccc_table );
                
              %% Round-Robin scheduling with Max-C  
                % not a proper scheduling algorithm (so won't be using it
                % anytime soon...)
                %[ ccc_output(t, rb, :), ccc_output_jd(t, rb, :), ~, ~, ~ ] = rr_max_c( num_users, combination_table(current_comb,:), all_signal_power(:, :, rb), all_signal_power_outer(:, :, :, rb), noise_power, ccc_table );
                
              %% Round-Robin schedulig with Max-C/I
                [ ccc_output_ci(t, rb, :), ccc_output_ci_jd(t, rb, :), connection(t, rb, :), connection_jd(t, rb, :) ] = rr_max_ci( num_users, combination_table(current_comb,:), all_signal_power(:, :, rb), all_signal_power_outer(:, :, :, rb), noise_power, ccc_table );
                
              %% increment
                current_comb = current_comb + 1;
                if current_comb > tot_combinations
                    current_comb = 1;
                end
                
                current_user = current_user + 1;
                if current_user > num_users
                    current_user = 1;
                end
            end
        end
        toc
        
        throughput_one_user(drop, trial, :, :) = ccc_output_one_user;
%         throughput(drop, trial, :, :, :) = ccc_output;
%         throughput_jd(drop, trial, :, :, :) = ccc_output_jd;
        throughput_ci(drop, trial, :, :, :) = ccc_output_ci;
        throughput_ci_jd(drop, trial, :, :, :) = ccc_output_ci_jd;
        
    end
    
end

% output for now:
sum(sum(sum(sum(throughput_one_user)))) / num_drops / trial_per_drop
% sum(sum(sum(sum(sum(throughput))))) / num_drops / trial_per_drop
% sum(sum(sum(sum(sum(throughput_jd))))) / num_drops / trial_per_drop
sum(sum(sum(sum(sum(throughput_ci))))) / num_drops / trial_per_drop
sum(sum(sum(sum(sum(throughput_ci_jd))))) / num_drops / trial_per_drop



%% Metric:


