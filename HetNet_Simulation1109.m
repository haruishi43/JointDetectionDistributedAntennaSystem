clear;
load('CCCtable_2antenna');

%% Randomize Every Simulation %%
rng('shuffle');

%% シミュレーション条件 %%
NO_user = 3;                         % ユーザ数
% NO_Interference_cell = 19;           % 19 Hexagonal area
NO_cell = 7;                         % picoBS数

% Center_frequency = 2.0*10^9;         % 中心周波数
NO_SC_inRB = 12;                     % 1RBに含まれるsubcarrier数
NO_RB = 24;                          % 1OFDMシンボルに含まれるRBの数
NO_SC = NO_SC_inRB * NO_RB;          % 1OFDMシンボルに含まれるsubcarrierの数
Band_RB = 180*10^3;                  % 1RBに使用する周波数帯域
% Band_SC = Band_RB / NO_SC_inRB;      % 1Subcarrierに使用する周波数帯域
Band = Band_RB * NO_RB;              % 使用する周波数帯域
NO_time_trial = 3;                   % 時間の試行回数 (3)
Timing_interval = 60;                % チャネルを固定するインターバル                                                   
NO_drop_trial = 1200;                % ユーザドロップの試行回数 (1200)
TI = 60;                             % Time Interval
NO_path = 6;                         % Jake'sモデルにおけるパスの数
% Doppler = 5.55;                      % Jake'sモデルにおけるドップラーシフト値[Hz]
% Decay = 1;                           % Jake'sモデルにおけるパスごとの減衰量
Interval = 1 / Band;                 % Jake'sモデルにおけるサンプリングインターバル
% JI = 600;                            % JakeInterval...fade.mを一周させるために数サンプルごとに取り出すようにする
rms_delay_spread = 1.0 * 10^(-6);    % 遅延スプレッド値
IS_distance = 500;                   % Inter-site distance
Shadowing_ave = 0;                   % Shadowingの平均値
Shadowing_var_Macro = 8;             % MacroBS Shadowingの分散値（？）
Shadowing_var_Pico = 10;             % PicoBS Shadowingの分散値（？）
% Correlation_jake_coefficient = 0.44; % 相関係数0.5に対応してjakeに入れる係数
Receiver_noise_density = -174;
NO_EIRP_base = 1;                    % MacroBS放射電力制御係数
% NO_EIRP_picobase = 1;                % PicoBS放射電力制御係数(今は制御していない)


%% Saving to Folder %%

% remove a couple of data in the beginning 
remove_beginning = 20;
% we will use the actual_interval
actual_interval = Timing_interval - remove_beginning;
% trial count
trial_count = 1;

% for each NO_time_trial we save, 
% input data:
channel_response = zeros(NO_RB, NO_user, NO_cell);
past_throughput = zeros(NO_user, 1);
% label data:
combination = zeros(NO_RB, 1);
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

SINR_SC_user = zeros(NO_user,NO_SC,(NO_cell+1)^NO_user );
SINR_SC_user_floor = zeros(NO_user,NO_SC,(NO_cell+1)^NO_user);
EIRP_base = zeros(1,NO_EIRP_base);
SINR_RB_user = zeros(NO_user,NO_RB,NO_EIRP_base);
SINR_RB_user_floor = zeros (NO_user,NO_RB,NO_EIRP_base);

channel_responce_freq = zeros(NO_user,NO_cell,NO_SC); %
channel_responce_time_rayleigh = zeros(NO_user,NO_cell,NO_SC); %

Coordinates = zeros(1,NO_user);                                                      % ユーザの座標を格納
% Coordinates_Interference = zeros(1,NO_Interference_cell);                            % 干渉基地局(Macro)の座標を格納                     

Capacity_byuser_macro_Conv = zeros(NO_user,NO_EIRP_base,NO_drop_trial);                                   
Delay_profile = zeros(1,NO_path); % 遅延

Distance_fromBS = zeros(NO_drop_trial,NO_cell,NO_user);
Distance_fromBS_pre = zeros(NO_drop_trial,NO_cell,NO_user);
PLR_fromBS = zeros(NO_drop_trial,NO_cell,NO_user);

Signal_power_fromBS_user = zeros(NO_user,NO_SC,NO_EIRP_base,(NO_cell+1)^NO_user);
Signal_power_fromBS_Interference_user = zeros(NO_user,NO_SC,NO_EIRP_base,(NO_cell+1)^NO_user);

%% BSの座標決定 %%
Coordinates_antenna = zeros(NO_cell);
Coordinates_antenna(1) = 0;

%% Coordinate Testing %%
Preset_Coordinates = [7, 1, 2]; % if needed

for a = 2:7
    Coordinates_antenna(a) = IS_distance * cos(a * pi/3 - pi/6) + 1i * IS_distance * sin(a * pi/3 - pi/6);
end

% なんで1-6だけ使っているのか？そもそもCoordinates_Interference使われていない
% for BS_index = 1:6
%     Coordinates_Interference(BS_index) = 2*IS_distance * cos(BS_index * pi/3 - pi/6) + 2i * IS_distance * sin(BS_index * pi/3 - pi/6)+IS_distance * cos(BS_index * pi/3 - pi/6 + pi/3) + 1i * IS_distance * sin(BS_index * pi/3 - pi/6 + pi/3);
%     for a = 1:6
%         Coordinates_Interference(BS_index * 6 + a) = Coordinates_Interference(BS_index) + IS_distance * cos(a * pi/3 - pi/6) + 1i * IS_distance * sin(a * pi/3 - pi/6);
%     end
% end

%% 遅延プロファイル作成 %%
for path_index = 1:NO_path
    Delay_profile(path_index) = exp( - (path_index-1) / (rms_delay_spread/Interval) );
end
Delay_profile = Delay_profile / sum(Delay_profile);

%% ユーザ配置を変えてシミュレーション %%    
for Drop_index = 1:NO_drop_trial
    tic
    disp(Drop_index);

    %% ユーザ配置の決定 %%
    for user_index = 1:NO_user
        Coordinates(user_index) = 0;        % initialization
        while Coordinates(user_index) == 0
            dx = (rand-0.5) * 2 * IS_distance / sqrt(3);
            dy = (rand-0.5) * 2 * IS_distance / sqrt(3);
            if (abs(dx) - IS_distance / 2 / sqrt(3)) * sqrt(3) > IS_distance / 2 - abs(dy) || abs(dy) > IS_distance / 2 ... %セル半径289mの六角形
                || abs(dx+dy*1i) < 10 %マクロ基地局とユーザの最小距離   
                
                Coordinates(user_index) = 0;
            else
                user_cell = randi(7);
                %user_cell = Preset_Coordinates(user_index);
%                 if user_cell == 1
%                     Coordinates(user_index) = dx + dy*1i;
%                 else
%                     % Coordinates(user_index) = Coordinates_antenna(user_cell - 1) + dx + dy*1i;
%                     % this part was wrong
%                     % change to: Coordinates(user_index) = Coordinates_antenna(user_cell) + dx + dy*1i;
%                 end
                Coordinates(user_index) = Coordinates_antenna(user_cell) + dx + dy*1i;
            end
        end
    end

    for cell_index = 1:NO_cell
        Distance_fromBS_pre(Drop_index,cell_index,:) = abs(Coordinates - repmat(Coordinates_antenna(cell_index),1,NO_user));     %　マクロ基地局からの直線距離
        Distance_fromBS(Drop_index,cell_index,:) = abs(sqrt(Distance_fromBS_pre(Drop_index,cell_index,:).^2 + 8.5^2));           %　マクロ基地局からの3D距離
        PLR_fromBS(Drop_index,cell_index,:) = 140.7 + 36.7 * log10(Distance_fromBS(Drop_index,cell_index,:)*0.001);  
    end

    %% アンテナパターン %%
    % ここはどこで使える？
    Rank_distance = zeros(1,NO_user);
    Coordinates_pre = Coordinates;
    for Divide_index = 1:NO_user/2
        for user_index = 1:NO_user
            if abs(Coordinates_pre(user_index)) == max(abs(Coordinates_pre))
                Rank_distance(user_index) = -1;
                Coordinates_pre(user_index) = 0;
                break
            end
        end
    end
    pow_amp = [10^(-0.30) 10^(-0.00) 10^(-0.20)  10^(-0.6)  10^(-0.8)  10^(-1.0)];
    tot_pow = sum(pow_amp);
     
    
    %% 時間ごとに容量計算 %%
    for Trial_time_index = 1:NO_time_trial
        Capacity_band_ave_macro_Conv = zeros(1,NO_user);    % Band毎の平均Capacityを格納
    
       %% Rayleigh Fading付加 %%
        GI = 32;% not used    
        for user_index = 1:NO_user
            for cell_index = 1:NO_cell
                channel_responce_time_rayleigh(user_index,cell_index,1:length(pow_amp)) = (1/sqrt(2).*(randn(1,length(pow_amp))+1j*randn(1,length(pow_amp)))) .* sqrt(Delay_profile);
                channel_responce_time_rayleigh(user_index,cell_index,7:NO_SC) = zeros(1,1,NO_SC-NO_path);
                channel_responce_freq(user_index,cell_index,:) = fft(channel_responce_time_rayleigh(user_index,cell_index,:));
            end
        end
        
        %% Average （→ チャネル応答）%%
        for user_index = 1:NO_user
            for cell_index = 1:NO_cell
               for RB_index = 1:NO_RB
                    channel_response(RB_index, user_index, cell_index) = abs(mean(channel_responce_freq(user_index, cell_index, NO_SC_inRB*(RB_index-1)+1:NO_SC_inRB*RB_index)));
               end
            end
        end
        
       %% シャドウイング計算 %%
       % 周辺BSからの干渉は考えられていない
%         Shadowing_correlation_matrix=ones(42,42);
%         Shadowing_correlation_matrix=Shadowing_correlation_matrix*0.5;
%         Shadowing=zeros(42,1);
%         for BS_index = 1:42
%             Shadowing(BS_index)=Shadowing_var_Pico.*randn(1,1)+Shadowing_ave; % Shadowing_ave is always 0? % Why Shadowing_var_pico
%             Shadowing_correlation_matrix(BS_index,BS_index)=1;
%         end
%         Shadowing_correlation_matrix=sqrtm(Shadowing_correlation_matrix);
%         Shadowing=Shadowing_correlation_matrix*Shadowing;
%         Shadowing=10.^((Shadowing)/10);
        
        % 本来ユーザごとに異なるはずなのでコメントアウトした
%         Shadowing_macro = sqrt(Shadowing_var_Macro).*randn(1,1) + Shadowing_ave; % log normalization (どこをとるかはrand)
%         Shadowing_macro = 10.^((Shadowing_macro)/10);
        
        %% Antenna Pairing %%
        a = zeros(1,user_index);
        notpair = zeros(1,(NO_cell + 1)^NO_user);
        bb=1;
        for user_antenna_pair = 1:(NO_cell+1)^NO_user
            user_antenna_pair_shift = user_antenna_pair-1;
            for user_index = 1:NO_user
                a(user_index) = fix(user_antenna_pair_shift/(NO_cell+1)^(NO_user-user_index))+1;
                user_antenna_pair_shift = rem(user_antenna_pair_shift,(NO_cell+1)^(NO_user-user_index));
            end
            if a(1) ~= 8
                if a(1) == a(2)
                    notpair(bb) = user_antenna_pair;
                    bb = bb+1;
                elseif a(1) == a(3)
                    notpair(bb) = user_antenna_pair;
                    bb = bb+1;
                end
            end
            if a(2) ~= 8
                if a(2) == a(3)
                    notpair(bb) = user_antenna_pair;
                    bb = bb+1;
                end
            end
        end
            
        
        %% SINR計算 %%
        Noise_power = 1;   
        EIRP_index = 1;
        EIRP_base(EIRP_index) = 0 + 30 - (Receiver_noise_density + 10*log10(Band));       
        for user_antenna_pair = 1:(NO_cell+1)^NO_user 
            nono = find(user_antenna_pair == notpair(:));
            if isempty(nono) == 1
                user_antenna_pair_shift = user_antenna_pair - 1; % なんで1引く？
                for user_index = 1:NO_user     
                    cell_index_select = fix(user_antenna_pair_shift/(NO_cell+1)^(NO_user-user_index))+1;
                    user_antenna_pair_shift = rem(user_antenna_pair_shift,(NO_cell+1)^(NO_user-user_index));

                    if cell_index_select == 8
                        Signal_power_fromBS_user(user_index,:,EIRP_index) = 0;
                    else
                        % シャドーイングの値をsqrt(Shadowing_var_Macro).*randn(1,1)とし，毎ユーザで値を変える
                        % 基地局・ユーザの干渉
                        Signal_power_fromBS_user(user_index,:,EIRP_index,user_antenna_pair) = sqrt(Shadowing_var_Macro).*randn(1,1) * 10.^((EIRP_base(EIRP_index)  - repmat(PLR_fromBS(Drop_index,cell_index_select,user_index)',NO_SC,1))/10) .* (abs(squeeze(channel_responce_freq(user_index,cell_index_select,:))).^2);    %マクロBSからの受信電力 
                        for user_index_int = 1:NO_user
                            if user_index_int ~= user_index
                                % シャドーイングの値をsqrt(Shadowing_var_Macro).*randn(1,1)とし，毎ユーザで値を変える
                                % ユーザ同士の干渉
                                Signal_power_fromBS_Interference_user(user_index_int,:,EIRP_index,user_antenna_pair) = squeeze(Signal_power_fromBS_Interference_user(user_index_int,:,EIRP_index,user_antenna_pair)).' + abs(sqrt(Shadowing_var_Macro).*randn(1,1) * 10.^((EIRP_base(EIRP_index) - repmat(PLR_fromBS(Drop_index,cell_index_select,user_index_int)',NO_SC,1))/10) .* (abs(squeeze(channel_responce_freq(user_index_int,cell_index_select,:))).^2));    %マクロBSからの受信電力
                            end
                        end
                    end
                end
                for user_index = 1:NO_user
                    SINR_SC_user(user_index,:,user_antenna_pair) = 10*log10(abs(Signal_power_fromBS_user(user_index,:,EIRP_index,user_antenna_pair)) ./ (Noise_power + Signal_power_fromBS_Interference_user(user_index,:,EIRP_index,user_antenna_pair)));
                end
            else
                SINR_SC_user(:,:,user_antenna_pair) = -inf;
            end
        end
        
        EIRP_index = 1;
        for user_antenna_pair = 1:(NO_cell+1)^NO_user
            for user_index = 1:NO_user
               for RB_index = 1:NO_RB
                   % average
                    SINR_RB_user(user_index,RB_index,user_antenna_pair) = mean(SINR_SC_user(user_index,NO_SC_inRB*(RB_index-1)+1:NO_SC_inRB*RB_index,user_antenna_pair));
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
        for user_antenna_pair = 1:(NO_cell+1)^NO_user
            for user_index = 1:NO_user
                for SC_index = 1:NO_SC
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
        
        Max_PFmetric_RB_Conv = zeros(Timing_interval,NO_RB);
        Capacity_Analyze_band_ave_macro_Conv = zeros(Timing_interval,NO_user); % ？
        select_user_antenna_pair = zeros(1,NO_RB);
        
        for Timing_interval_index = 1:Timing_interval
            %% PFmetric計算
            Max_CC_modulation = zeros(NO_user,NO_RB,(NO_cell+1)^NO_user);
            modulation_index = zeros(NO_user,NO_RB,(NO_cell+1)^NO_user);
            select_usermod = zeros(NO_RB,NO_user);
            for RB_index = 1:NO_RB
                EIRP_index = 1;
                PFmetric_RB_Conv = zeros(1,(NO_cell+1)^NO_user);
                for user_antenna_pair = 1:(NO_cell+1)^NO_user
                    if Capacity_band_ave_macro_Conv(user_index) == 0
                        Conv_Capacity_band_ave_pre  = 1;
                    else
                        Conv_Capacity_band_ave_pre = Capacity_band_ave_macro_Conv(user_index);
                    end
                    for user_index = 1:NO_user
                        if SINR_RB_user_floor(user_index,RB_index,user_antenna_pair) ~= -inf
                            [Max_CC_modulation(user_index,RB_index,user_antenna_pair),modulation_index(user_index,RB_index,user_antenna_pair)] = max(squeeze(CCCtable_conv_SINRp_alphap_QAMq_QAMp(SINR_RB_user_floor(user_index,RB_index,user_antenna_pair)+11,1,1,:)));
                            PFmetric_RB_Conv(user_antenna_pair) = PFmetric_RB_Conv(user_antenna_pair) + Max_CC_modulation(user_index,RB_index,user_antenna_pair)*NO_SC*Timing_interval*NO_time_trial / Conv_Capacity_band_ave_pre; %%Conv_Capacity_band_ave_preが今までの平均通信路容量
                        end
                    end
                end
                if max(PFmetric_RB_Conv) >= Max_PFmetric_RB_Conv(Timing_interval_index,RB_index)
                    Max_PFmetric_RB_Conv(Timing_interval_index,RB_index) = max(PFmetric_RB_Conv);
                    index = find(PFmetric_RB_Conv == Max_PFmetric_RB_Conv(Timing_interval_index,RB_index));
                    % what I want:
                    select_user_antenna_pair(RB_index) = index(1,1);
                    for user_index = 1:NO_user
                        select_usermod(RB_index,user_index) = modulation_index(user_index,RB_index,select_user_antenna_pair(RB_index));
                    end
                end
            end
            
            
            %% 容量計算 (従来方式)
            Capacity_trial_realnear_prop = 0;
            Capacity_trial_realfar_prop = 0;
            capacity_macro_Conv_SC = zeros(NO_user,NO_SC);
            Capacity_byuser_cur_macro_Conv = zeros(1,NO_user);

            for SC_index = 1:NO_SC
                RB_index = floor((SC_index-1)/NO_SC_inRB)+1;
                for user_index = 1:NO_user
                    SINR_user = SINR_SC_user_floor(user_index,SC_index,select_user_antenna_pair(RB_index))+11;
                    if SINR_user ~= -inf
                        farmod_index = select_usermod(RB_index,user_index);
                        capacity_macro_Conv_SC(user_index,SC_index) = CCCtable_conv_SINRp_alphap_QAMq_QAMp(SINR_user,1,1,farmod_index);     %所望信号容量
                        
                        %% add capacity ave.
                        Capacity_byuser_cur_macro_Conv(user_index) = Capacity_byuser_cur_macro_Conv(user_index) + capacity_macro_Conv_SC(user_index,SC_index) ;
                        
                        if Timing_interval_index > remove_beginning
                            Capacity_byuser_macro_Conv(user_index,1,Drop_index) = Capacity_byuser_macro_Conv(user_index,1,Drop_index) + capacity_macro_Conv_SC(user_index,SC_index) ;
                        end
                    end
                end
            end
            
            %% saving data
            % at 60
            if Timing_interval_index == 60
                % Get combination (base 10 to base 8)
                for rb_index= 1:NO_RB
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
            Capacity_band_ave_macro_Conv = (1-1/TI)*Capacity_band_ave_macro_Conv + 1/TI*Capacity_byuser_cur_macro_Conv;
            % all: 使われていない
            Capacity_Analyze_band_ave_macro_Conv(Timing_interval_index,:) = Capacity_band_ave_macro_Conv;
            
            
        end
        % end of timing_interval
        
        %% saving data
        %channel_response
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
    disp("end of trials");
    toc
end