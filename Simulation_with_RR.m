clear;
load('CCCtable_2antenna');

%% Randomiz:
rng('shuffle');

%% Model parameters:
num_users = 3;                      % # of users
num_cell = 7;                       % # of cell

num_rb = 24;                        % # of resource blocks in 1 OFDM symbol
num_sc_in_rb = 12;                  % # of subcarriers in resource blocks
num_sc = num_rb * num_sc_in_rb;     % # of total subcarriers

band_per_rb = 180*10^3;             % frequency band range for each rb (Hz)
band = band_per_rb * num_rb;        % total frequency band

rnd = -174;                          % receiver noise density
noise_power = 1;
eirp = 0 + 30 - ( rnd + 10*log10(band) );

%% Simulation parameters:
num_drops = 30;
time_interval = 60;
trial_per_drop = 1;

%% Initializing variables:
plr_from_bs_all = zeros(num_drops,num_cell,num_users);              % packet loss ratio 
channel_response_freq = zeros(num_users, num_cell, num_sc);         % channel response frequency
[pairs, impossible_pairs] = create_pairs( num_users, num_cell );                        % antenna pair list
num_pairs = numel(pairs);   % # of pairs

%% シミュレーション条件 %% old

Shadowing_ave = 0;                   % Shadowingの平均値
Shadowing_var_Macro = 8;             % MacroBS Shadowingの分散値（？）
Shadowing_var_Pico = 10;             % PicoBS Shadowingの分散値（？）

NO_EIRP_base = 1;                    % MacroBS放射電力制御係数


%% Saving to Folder %%

% remove a couple of data in the beginning 
remove_beginning = 20;
% we will use the actual_interval
actual_interval = time_interval - remove_beginning;
% trial count
trial_count = 1;

% for each NO_time_trial we save, 
% input data:
channel_response = zeros(num_rb, num_users, num_cell);
past_throughput = zeros(num_users, 1);
% label data:
combination = zeros(num_rb, 1);
% other data:
sumrate = 0;
% which means we get actual_interval as the number of data points per each
% NO_time_trail

% for saving, we divide the data for each NO_time_trail

% create the main folder for saving:
data_folder_name = datestr(datetime('now', 'TimeZone','local','Format','y-MM-dd_HH:mm:ss'), 'yyyy-mm-dd_HH-MM-SS');
[status, msg] = mkdir(data_folder_name);
if status ~= 1
    disp(msg);
end
% change directory:
% home: starting directory
%home = cd(data_folder_name); 
home = pwd;
data_root = fullfile(home, data_folder_name);

%% More Variables 

EIRP_base = zeros(1,NO_EIRP_base);

Capacity_byuser_macro_Conv = zeros(num_users,NO_EIRP_base,num_drops);                                   

signal_pow_from_bs = zeros( num_users, num_sc, num_pairs );
signal_pow_from_users = zeros( num_users, num_sc, num_pairs );
sinr_user_sc = zeros(num_users, num_sc, num_pairs);
sinr_user_sc_floor = zeros(num_users, num_sc, num_pairs);
sinr_user_rb = zeros(num_users, num_rb, num_pairs);
sinr_user_rb_floor = zeros(num_users, num_rb, num_pairs);

%% Create coordinates for each BS:
preset_coordinates = [1 2 7]; % For Coordinate Testing
antenna_coordinates = create_bs_coordinate();

%% Simulation loop (change user placement):   
for drop = 1:num_drops
    tic
    
    %% Create Coordinates for each user:
    user_coordinates = create_user_coordinates( antenna_coordinates );
    
    %% Calculate Packet Loss Ratio
    plr_from_bs_all(drop, :, :) = create_plr_from_bs( antenna_coordinates, user_coordinates );
    
    %% Simulation loop (trial per drop):
    for trial = 1:trial_per_drop
        
        %% Calculate Rayleigh Fading:
        channel_response_freq = add_rayleigh_fading( num_users, num_cell );
        
        %% Average to create channel response (for saving purpose):
        for user = 1:num_users
            for cell = 1:num_cell
               for rb = 1:num_rb
                    channel_response(rb, user, cell) = abs(mean(channel_response_freq(user, cell, num_sc_in_rb * (rb-1) + 1:num_sc_in_rb * rb)));
               end
            end
        end 
        
        %% Calculate SINR:
        for pair_index = 1:num_pairs
            pair_shift = pairs(pair_index) - 1;
            
            for i = 1:num_users
                selected_cell = fix(pair_shift / (num_cell + 1)^(num_users - i)) + 1;
                pair_shift = rem( pair_shift, (num_cell + 1)^(num_users - i) );
                
                if selected_cell == 8
                    % zero if rest
                    signal_pow_from_bs(i, :, pair_index) = 0;
                else
                    % シャドーイングの値をsqrt(Shadowing_var_Macro).*randn(1,1)とし，毎ユーザで値を変える
                    % 基地局・ユーザの干渉
                    signal_pow_from_bs(i, :, pair_index) = sqrt(Shadowing_var_Macro).*randn(1,1) * 10.^((eirp - repmat(plr_from_bs_all(drop, selected_cell, i)', num_sc,1))/10) .* (abs(squeeze(channel_response_freq(i,selected_cell,:))).^2);
                    
                    for j = 1:num_users
                        if i ~= j
                            % シャドーイングの値をsqrt(Shadowing_var_Macro).*randn(1,1)とし，毎ユーザで値を変える
                            % ユーザ同士の干渉
                            signal_pow_from_users(j, :, pair_index) = squeeze(signal_pow_from_users(j,:,pair_index)).' + abs(sqrt(Shadowing_var_Macro).*randn(1,1) * 10.^((eirp - repmat(plr_from_bs_all(drop,selected_cell,j)',num_sc,1))/10) .* (abs(squeeze(channel_response_freq(j,selected_cell,:))).^2));
                        end
                    end
                end
            end
            
            for i = 1:num_users
                
                sinr_user_sc(i, :, pair_index) = 10*log10(abs(signal_pow_from_bs(i,:,pair_index)) ./ (noise_power + signal_pow_from_users(i,:,pair_index)) );
            
                % Fix SINR between -10 and 30:
                % For Subcarriers
                for sc = 1:num_sc
                    sinr_user_sc_floor(i,sc,pair_index) = floor(sinr_user_sc(i,sc,pair_index));
                    if sinr_user_sc_floor(i,sc,pair_index) <= -10
                        sinr_user_sc_floor(i,sc,pair_index) = -10;
                    elseif sinr_user_sc_floor(i,sc,pair_index) >= 30
                        sinr_user_sc_floor(i,sc,pair_index) = 30;
                    end
                end
                
                % For resource blocks
                for rb = 1:num_rb
                    sinr_user_rb(i,rb,pair_index) = mean(sinr_user_sc(i,num_sc_in_rb*(rb - 1)+1:num_sc_in_rb*rb,pair_index));
                    sinr_user_rb_floor(i,rb,pair_index) = floor(sinr_user_rb(i,rb,pair_index));
                    if sinr_user_rb_floor(i,rb,pair_index) <= -10
                        sinr_user_rb_floor(i,rb,pair_index) = -10;
                    elseif sinr_user_rb_floor(i,rb,pair_index) >= 30
                        sinr_user_rb_floor(i,rb,pair_index) = 30;
                    end
                end
                
            end
        end
        
        %% Simulation loop (for time interval):
        Capacity_band_ave_macro_Conv = zeros(1, num_users);    % Band毎の平均Capacityを格納
        Max_PFmetric_RB_Conv = zeros(time_interval,num_rb);
        selected_pair_for_rb = zeros(1, num_rb);
        
        for t = 1:time_interval
            %% Calculate PF metric:
            Max_CC_modulation = zeros(num_users,num_rb,num_pairs);
            modulation_index = zeros(num_users,num_rb,num_pairs);
            select_usermod = zeros(num_rb,num_users);
            for rb = 1:num_rb
                EIRP_index = 1;
                PFmetric_RB_Conv = zeros(1,num_pairs);
                for pair_index = 1:num_pairs
                    if Capacity_band_ave_macro_Conv(i) == 0
                        Conv_Capacity_band_ave_pre  = 1;
                    else
                        Conv_Capacity_band_ave_pre = Capacity_band_ave_macro_Conv(i);
                    end
                    for i = 1:num_users
                        %%% why 11?
                        [Max_CC_modulation(i,rb,pair_index),modulation_index(i,rb,pair_index)] = max(squeeze(CCCtable_conv_SINRp_alphap_QAMq_QAMp(sinr_user_rb_floor(i,rb,pair_index)+11,1,1,:))); 
                        PFmetric_RB_Conv(pair_index) = PFmetric_RB_Conv(pair_index) + Max_CC_modulation(i,rb,pair_index)*num_sc*time_interval*trial_per_drop / Conv_Capacity_band_ave_pre; %%Conv_Capacity_band_ave_preが今までの平均通信路容量
                    end
                end
                if max(PFmetric_RB_Conv) >= Max_PFmetric_RB_Conv(t,rb)
                    Max_PFmetric_RB_Conv(t,rb) = max(PFmetric_RB_Conv);
                    index = find(PFmetric_RB_Conv == Max_PFmetric_RB_Conv(t,rb));
                    selected_pair_for_rb(rb) = index(1,1);
                    for i = 1:num_users
                        select_usermod(rb,i) = modulation_index(i,rb,selected_pair_for_rb(rb));
                    end
                end
            end
            
            %% 容量計算 (従来方式)
            capacity_macro_Conv_SC = zeros(num_users,num_sc);
            Capacity_byuser_cur_macro_Conv = zeros(1,num_users);

            for sc = 1:num_sc
                rb = floor((sc-1)/num_sc_in_rb)+1;
                for i = 1:num_users
                    SINR_user = sinr_user_sc_floor(i,sc,selected_pair_for_rb(rb))+11;
                    if SINR_user ~= -inf
                        farmod_index = select_usermod(rb,i);
                        capacity_macro_Conv_SC(i,sc) = CCCtable_conv_SINRp_alphap_QAMq_QAMp(SINR_user,1,1,farmod_index);     %所望信号容量
                        
                        %% add capacity ave.
                        Capacity_byuser_cur_macro_Conv(i) = Capacity_byuser_cur_macro_Conv(i) + capacity_macro_Conv_SC(i,sc) ;
                        
                        if t > remove_beginning
                            Capacity_byuser_macro_Conv(i,1,drop) = Capacity_byuser_macro_Conv(i,1,drop) + capacity_macro_Conv_SC(i,sc) ;
                        end
                    end
                end
            end
            
            %% saving data
            % at 60
            if t == 60
                % Get combination (base 10 to base 8)
                for rb_index= 1:num_rb
                    a1 = floor((selected_pair_for_rb(1, rb_index)-1)/64) + 1;
                    tmp = mod((selected_pair_for_rb(1, rb_index)-1), 64);
                    a2 = floor(tmp/8) + 1;
                    tmp = mod(tmp, 8);
                    a3 = floor(tmp/1) + 1;
                    combination(rb_index, 1) = 100*a1 + 10*a2 + a3;
                end
                
                % Get past throughputs
                past_throughput(:, 1) = Capacity_band_ave_macro_Conv.';
                
                % Get sumrate of the network
                sumrate = sum(Capacity_byuser_cur_macro_Conv); 
            end
            
            % one:
            Capacity_band_ave_macro_Conv = (1-1/time_interval)*Capacity_band_ave_macro_Conv + 1/time_interval*Capacity_byuser_cur_macro_Conv;

        end
        % end of timing_interval
        
        %% saving data
        dir_trial_name = int2str(trial_count);
        basename = fullfile(data_root, dir_trial_name);
        [status, msg] = mkdir(basename);
        if status ~= 1
            disp(msg); % if any error
        end
        cd(basename);
        
        save('channel_response.mat', 'channel_response');
        save('combination.mat', 'combination');
        save('past_throughput.mat', 'past_throughput');
        save('sumrate.mat', 'sumrate');
        
        % clean up:
        trial_count = trial_count + 1;
        cd(home)
        
    end 
    % end of trials
    %disp("end of trials");
    toc
end