clear;
ccc_table = load('CCCtable_2antenna', ...
                 'CCCtable_conv_SINRp_alphap_QAMq_QAMp', ...   % no joint ml detection
                 'CCCtable_prop_SINRp_alphap_QAMq_QAMp');      % joint ml detection

%% Randomize:
rng('Shuffle');

%% Model parameters:
users = [5 10 15 20];
user_trial = numel(users);
num_cell = 7;                       % # of cell
num_outer_macro = 6;
preset_coordinates = [1 3 5 6 7];   % For Coordinate Testing (has to change when num_users change)

num_rb = 24;                        % # of resource blocks in 1 OFDM symbol
num_sc_in_rb = 12;                  % # of subcarriers in resource blocks
num_sc = num_rb * num_sc_in_rb;     % # of total subcarriers

distance = 25;

band_per_rb = 180*10^3;             % frequency band range for each rb (Hz)
band = band_per_rb * num_rb;        % total frequency band

shadowing_ave = 0;
shadowing_var = 8;
rnd = -174;                         % Reciever Noise Density
noise_power = rnd + 10*log10( band );
eirp = 0 + 30;

% Scheduling parameters
num_select = 2;                     % # of user selected for each combination

%% Simulation parameters:
num_drops = 100;
trial_per_drop = 2;
time_interval = 30;

%% Saving variables:
throughput_one_user = zeros(user_trial, num_drops, trial_per_drop, time_interval, num_rb);
% throughput = zeros(user_trial, num_drops, trial_per_drop, time_interval, num_rb, num_select);
% throughput_jd = zeros(user_trial, num_drops, trial_per_drop, time_interval, num_rb, num_select);
throughput_ci = zeros(user_trial, num_drops, trial_per_drop, time_interval, num_rb, num_select);
throughput_ci_jd = zeros(user_trial, num_drops, trial_per_drop, time_interval, num_rb, num_select);

% 1: 1 user scheduling
% 2: max-ci without jd
% 3: max-ci with jd
ave_throughput_num_users = zeros(user_trial, 3);

for user = 1:user_trial
    user
    %% Model parameters:
    num_users = users(user);                   % # of users
    
    % Scheduling Combinations
    [combination_table, tot_combinations] = create_combination_table( num_users, num_select );

    %% Initializing variables:
    plr_from_bs_all = zeros(num_drops, num_users, num_cell);              % propagation loss ratio
    plr_from_outer_cell = zeros(num_drops, num_outer_macro, num_users, num_cell);          % propagation loss ratio from outer cell

    channel_response_freq = zeros(num_users, num_cell, num_sc);           % channel response frequency
    channel_response = zeros(num_users, num_cell, num_rb);

    %% Create coordinates for each BS:
    antenna_coordinates = create_bs_coordinate( distance );

    %% Create outer cell coordinates:
    outer_cell_coordinates = create_outer_cell_coordinates( distance );

    %% Simulation loop (change user placement): 
    tic
    for drop = 1:num_drops

        %% Create Coordinates for each user:
        user_coordinates = create_user_coordinates( antenna_coordinates, num_users );

        %% Calculate Propagation Loss 
        plr_from_bs_all(drop, :, :) = create_plr_from_bs( antenna_coordinates, user_coordinates );
        plr_from_outer_cell(drop, :, :, :) = create_plr_from_outer_cell( outer_cell_coordinates, user_coordinates );

        %% Simulation loop (trial):
        for trial = 1:trial_per_drop

            %% Calculate Rayleigh Fading:
            channel_response_freq = add_rician_fading( num_users, num_cell );

            %% Average to create channel response for each RB:
            all_signal_power = zeros(num_users, num_cell, num_rb);
            for u = 1:num_users
                for cell = 1:num_cell

                    const = 10.^(( eirp  - plr_from_bs_all(drop, u, cell) ) / 10);

                    for rb = 1:num_rb

                        % channel response (average of all subcarriers in a
                        % resource block) 
                        channel_response(u, cell, rb) = mean( channel_response_freq( u, cell, num_sc_in_rb * (rb-1) + 1:num_sc_in_rb * rb ) );
                        % this might be better:
                        %channel_response(user, cell, rb) = mean( abs( channel_response_freq( user, cell, num_sc_in_rb * (rb-1) + 1:num_sc_in_rb * rb ) ).^2 );

                        % signal in real number domain
                        all_signal_power(u, cell, rb) = 10^( sqrt(shadowing_var)*randn(1,1) / 10 ) * const * ( abs( channel_response(u, cell, rb) ).^2 );

                    end
                end
            end

            %% Signal power from outer cell:
            all_signal_power_outer = zeros(num_outer_macro, num_users, num_cell, num_rb);
            for macro = 1:num_outer_macro
                channel_response_macro = zeros(num_users, num_cell, num_rb);
                channel_response_macro_freq = add_rayleigh_fading( num_users, num_cell );

                for u = 1:num_users
                    for cell = 1:num_cell

                        const = 10.^(( eirp  - plr_from_outer_cell(drop, macro, u, cell) ) / 10);

                        for rb = 1:num_rb

                            channel_response_macro(u, cell, rb) = mean( channel_response_macro_freq( u, cell, num_sc_in_rb * (rb-1) + 1:num_sc_in_rb * rb ) );
                            % signal in real number domain
                            all_signal_power_outer(macro, u, cell, rb) = const * ( abs( channel_response_macro(u, cell, rb) ).^2 );

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

            throughput_one_user(user, drop, trial, :, :) = ccc_output_one_user;
    %         throughput(user, drop, trial, :, :, :) = ccc_output;
    %         throughput_jd(user, drop, trial, :, :, :) = ccc_output_jd;
            throughput_ci(user, drop, trial, :, :, :) = ccc_output_ci;
            throughput_ci_jd(user, drop, trial, :, :, :) = ccc_output_ci_jd;

        end
        
    end
    toc
    % output for now:
    ave_throughput_num_users(user, 1) = sum(sum(sum(sum(throughput_one_user(user, :, :, :, :))))) / num_drops / trial_per_drop / num_rb / time_interval * (7/0.5 * 1000);
    % sum(sum(sum(sum(sum(throughput))))) / num_drops / trial_per_drop
    % sum(sum(sum(sum(sum(throughput_jd))))) / num_drops / trial_per_drop
    ave_throughput_num_users(user, 2) = sum(sum(sum(sum(sum(throughput_ci(user, :, :, :, :, :)))))) / num_drops / trial_per_drop / num_rb / time_interval * (7/0.5 * 1000);
    ave_throughput_num_users(user, 3) = sum(sum(sum(sum(sum(throughput_ci_jd(user, :, :, :, :)))))) / num_drops / trial_per_drop / num_rb / time_interval * (7/0.5 * 1000);
    
    ave_throughput_num_users
    
end

%% Metric:

%% Plot Average Throughput vs. Number of Users:

figure(1)
plot(users,ave_throughput_num_users(:, 1),'-y','LineWidth',3)
xlim([5 20]);
hold on
plot(users,ave_throughput_num_users(:, 2),'-.m','LineWidth',3);
hold on
plot(users,ave_throughput_num_users(:, 3),'-.c','LineWidth',3);
hold on
grid on

legend('Single-MT','Multi-MT w/o Joint Detection','Multi-MT with Joint Detection','Location','NorthWest')
xlabel('Number of MTs (N_{MT})','FontName','Arial','FontSize',14)
ylabel('Average Throughput (bit / RB / sec)','FontName','Arial','FontSize',14)
set(gca,'FontName','Arial','FontSize',10)
hold off


%% CDF

