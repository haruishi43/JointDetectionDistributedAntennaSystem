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

%% Simulation parameters:
num_drops = 30;
time_interval = 60;
trial_per_drop = 1;

%% Initializing variables:
plr_from_bs_all = zeros(num_drops,num_cell,num_users);              % packet loss ratio 
channel_response_freq = zeros(num_users, num_cell, num_sc);         % channel response frequency
impossible_list = create_impossible_pairs( num_users, num_cell );   % impossible antenna pair list



%% シミュレーション条件 %% old

Shadowing_ave = 0;                   % Shadowingの平均値
Shadowing_var_Macro = 8;             % MacroBS Shadowingの分散値（？）
Shadowing_var_Pico = 10;             % PicoBS Shadowingの分散値（？）
Receiver_noise_density = -174;
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

SINR_SC_user = zeros(num_users,num_sc,(num_cell+1)^num_users );
SINR_SC_user_floor = zeros(num_users,num_sc,(num_cell+1)^num_users);
EIRP_base = zeros(1,NO_EIRP_base);
SINR_RB_user = zeros(num_users,num_rb,NO_EIRP_base);
SINR_RB_user_floor = zeros (num_users,num_rb,NO_EIRP_base);

Capacity_byuser_macro_Conv = zeros(num_users,NO_EIRP_base,num_drops);                                   



Signal_power_fromBS_user = zeros(num_users,num_sc,NO_EIRP_base,(num_cell+1)^num_users);
Signal_power_fromBS_Interference_user = zeros(num_users,num_sc,NO_EIRP_base,(num_cell+1)^num_users);



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
        tic
        Capacity_band_ave_macro_Conv = zeros(1, num_users);    % Band毎の平均Capacityを格納
        
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
        
        
        %% SINR計算 %%
        Noise_power = 1;   
        EIRP_index = 1;
        EIRP_base(EIRP_index) = 0 + 30 - (Receiver_noise_density + 10*log10(band));       
        for user_antenna_pair = 1:(num_cell+1)^num_users 
            impossible = find(user_antenna_pair == impossible_list(:));
            if isempty(impossible) == 1
                user_antenna_pair_shift = user_antenna_pair - 1; % なんで1引く？
                for user_index = 1:num_users     
                    cell_index_select = fix(user_antenna_pair_shift/(num_cell+1)^(num_users-user_index))+1;
                    user_antenna_pair_shift = rem(user_antenna_pair_shift,(num_cell+1)^(num_users-user_index));

                    if cell_index_select == 8
                        Signal_power_fromBS_user(user_index,:,EIRP_index) = 0;
                    else
                        % シャドーイングの値をsqrt(Shadowing_var_Macro).*randn(1,1)とし，毎ユーザで値を変える
                        % 基地局・ユーザの干渉
                        Signal_power_fromBS_user(user_index,:,EIRP_index,user_antenna_pair) = sqrt(Shadowing_var_Macro).*randn(1,1) * 10.^((EIRP_base(EIRP_index)  - repmat(plr_from_bs_all(drop,cell_index_select,user_index)',num_sc,1))/10) .* (abs(squeeze(channel_response_freq(user_index,cell_index_select,:))).^2);    %マクロBSからの受信電力 
                        for user_index_int = 1:num_users
                            if user_index_int ~= user_index
                                % シャドーイングの値をsqrt(Shadowing_var_Macro).*randn(1,1)とし，毎ユーザで値を変える
                                % ユーザ同士の干渉
                                Signal_power_fromBS_Interference_user(user_index_int,:,EIRP_index,user_antenna_pair) = squeeze(Signal_power_fromBS_Interference_user(user_index_int,:,EIRP_index,user_antenna_pair)).' + abs(sqrt(Shadowing_var_Macro).*randn(1,1) * 10.^((EIRP_base(EIRP_index) - repmat(plr_from_bs_all(drop,cell_index_select,user_index_int)',num_sc,1))/10) .* (abs(squeeze(channel_response_freq(user_index_int,cell_index_select,:))).^2));    %マクロBSからの受信電力
                            end
                        end
                    end
                end
                for user_index = 1:num_users
                    SINR_SC_user(user_index,:,user_antenna_pair) = 10*log10(abs(Signal_power_fromBS_user(user_index,:,EIRP_index,user_antenna_pair)) ./ (Noise_power + Signal_power_fromBS_Interference_user(user_index,:,EIRP_index,user_antenna_pair)));
                end
            else
                SINR_SC_user(:,:,user_antenna_pair) = -inf;
            end
        end
        
        EIRP_index = 1;
        for user_antenna_pair = 1:(num_cell+1)^num_users
            for user_index = 1:num_users
               for RB_index = 1:num_rb
                   % average
                    SINR_RB_user(user_index,RB_index,user_antenna_pair) = mean(SINR_SC_user(user_index,num_sc_in_rb*(RB_index-1)+1:num_sc_in_rb*RB_index,user_antenna_pair));
                    SINR_RB_user_floor(user_index,RB_index,user_antenna_pair) = floor(SINR_RB_user(user_index,RB_index,user_antenna_pair));
                    if SINR_RB_user_floor(user_index,RB_index,user_antenna_pair) == -inf
                        SINR_RB_user_floor(user_index,RB_index,user_antenna_pair) = -inf;
                    elseif SINR_RB_user_floor(user_index,RB_index,user_antenna_pair) <=-10
                        SINR_RB_user_floor(user_index,RB_index,user_antenna_pair) = -10;
                    elseif SINR_RB_user_floor(user_index,RB_index,user_antenna_pair) >= 30
                        SINR_RB_user_floor(user_index,RB_index,user_antenna_pair) =30;
                    end
                end
            end   
        end

        %% SINR値を変域内に修正 %%
        for user_antenna_pair = 1:(num_cell+1)^num_users
            for user_index = 1:num_users
                for SC_index = 1:num_sc
                    SINR_SC_user_floor(user_index,SC_index,user_antenna_pair) = floor(SINR_SC_user(user_index,SC_index,user_antenna_pair));
                    if SINR_SC_user_floor(user_index,SC_index,user_antenna_pair) == -inf
                        SINR_SC_user_floor(user_index,SC_index,user_antenna_pair) = -inf;
                    elseif SINR_SC_user_floor(user_index,SC_index,user_antenna_pair) <= -10
                        SINR_SC_user_floor(user_index,SC_index,user_antenna_pair) = -10;
                    elseif SINR_SC_user_floor(user_index,SC_index,user_antenna_pair) >= 30
                        SINR_SC_user_floor(user_index,SC_index,user_antenna_pair) = 30;
                    end
                end
            end
        end
        
        Max_PFmetric_RB_Conv = zeros(time_interval,num_rb);
        %Capacity_Analyze_band_ave_macro_Conv = zeros(Timing_interval,NO_user); % ？
        select_user_antenna_pair = zeros(1,num_rb);
        
        for Timing_interval_index = 1:time_interval
            %% PFmetric計算
            Max_CC_modulation = zeros(num_users,num_rb,(num_cell+1)^num_users);
            modulation_index = zeros(num_users,num_rb,(num_cell+1)^num_users);
            select_usermod = zeros(num_rb,num_users);
            for RB_index = 1:num_rb
                EIRP_index = 1;
                PFmetric_RB_Conv = zeros(1,(num_cell+1)^num_users);
                for user_antenna_pair = 1:(num_cell+1)^num_users
                    if Capacity_band_ave_macro_Conv(user_index) == 0
                        Conv_Capacity_band_ave_pre  = 1;
                    else
                        Conv_Capacity_band_ave_pre = Capacity_band_ave_macro_Conv(user_index);
                    end
                    for user_index = 1:num_users
                        if SINR_RB_user_floor(user_index,RB_index,user_antenna_pair) ~= -inf
                            [Max_CC_modulation(user_index,RB_index,user_antenna_pair),modulation_index(user_index,RB_index,user_antenna_pair)] = max(squeeze(CCCtable_conv_SINRp_alphap_QAMq_QAMp(SINR_RB_user_floor(user_index,RB_index,user_antenna_pair)+11,1,1,:)));
                            PFmetric_RB_Conv(user_antenna_pair) = PFmetric_RB_Conv(user_antenna_pair) + Max_CC_modulation(user_index,RB_index,user_antenna_pair)*num_sc*time_interval*trial_per_drop / Conv_Capacity_band_ave_pre; %%Conv_Capacity_band_ave_preが今までの平均通信路容量
                        end
                    end
                end
                if max(PFmetric_RB_Conv) >= Max_PFmetric_RB_Conv(Timing_interval_index,RB_index)
                    Max_PFmetric_RB_Conv(Timing_interval_index,RB_index) = max(PFmetric_RB_Conv);
                    index = find(PFmetric_RB_Conv == Max_PFmetric_RB_Conv(Timing_interval_index,RB_index));
                    % what I want:
                    select_user_antenna_pair(RB_index) = index(1,1);
                    for user_index = 1:num_users
                        select_usermod(RB_index,user_index) = modulation_index(user_index,RB_index,select_user_antenna_pair(RB_index));
                    end
                end
            end
            
            
            %% 容量計算 (従来方式)
            %Capacity_trial_realnear_prop = 0;
            %Capacity_trial_realfar_prop = 0;
            capacity_macro_Conv_SC = zeros(num_users,num_sc);
            Capacity_byuser_cur_macro_Conv = zeros(1,num_users);

            for SC_index = 1:num_sc
                RB_index = floor((SC_index-1)/num_sc_in_rb)+1;
                for user_index = 1:num_users
                    SINR_user = SINR_SC_user_floor(user_index,SC_index,select_user_antenna_pair(RB_index))+11;
                    if SINR_user ~= -inf
                        farmod_index = select_usermod(RB_index,user_index);
                        capacity_macro_Conv_SC(user_index,SC_index) = CCCtable_conv_SINRp_alphap_QAMq_QAMp(SINR_user,1,1,farmod_index);     %所望信号容量
                        
                        %% add capacity ave.
                        Capacity_byuser_cur_macro_Conv(user_index) = Capacity_byuser_cur_macro_Conv(user_index) + capacity_macro_Conv_SC(user_index,SC_index) ;
                        
                        if Timing_interval_index > remove_beginning
                            Capacity_byuser_macro_Conv(user_index,1,drop) = Capacity_byuser_macro_Conv(user_index,1,drop) + capacity_macro_Conv_SC(user_index,SC_index) ;
                        end
                    end
                end
            end
            
            %% saving data
            % at 60
            if Timing_interval_index == 60
                % Get combination (base 10 to base 8)
                for rb_index= 1:num_rb
                    a1 = floor((select_user_antenna_pair(1, rb_index)-1)/64) + 1;
                    tmp = mod((select_user_antenna_pair(1, rb_index)-1), 64);
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
            % all: 使われていない
            %Capacity_Analyze_band_ave_macro_Conv(Timing_interval_index,:) = Capacity_band_ave_macro_Conv;
            
            
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