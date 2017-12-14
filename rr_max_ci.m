function [ ccc_output, ccc_output_jd, connection, connection_jd ] = rr_max_ci( num_users, combination, current_signal, all_signal_power_outer, noise, ccc_table )
% Round-Robin scheduling with Max-C/I scheduling 

% update best throughput modulation

%% Outputs
ccc_output = zeros(1, 2);
ccc_output_jd = zeros(1, 2);
connection = 8 * ones(1, num_users);
connection_jd = 8 * ones(1, num_users);

%% RR and Max-C/I
u1 = combination(1);
u2 = combination(2);

% signal (real domain) table for each user
current_signal_u1 = current_signal(u1, :);
current_signal_u2 = current_signal(u2, :);

best_throughput = 0;
best_throughput_jd = 0;
best_mod = [ 0 0 ];
best_mod_jd = [ 0 0 ];
cell_pair = [ 0 0 ];
cell_pair_jd = [ 0 0 ];
best_pow = [ 0 0 ];
best_pow_jd = [ 0 0 ]; 
best_alpha = [ 0 0 ];
best_alpha_jd = [ 0 0 ];

for cell_u1 = 1:7
    % variables
    signal = zeros(2, 2);
    power = zeros(1, 2);
    
    signal(1, 1) = current_signal_u1(cell_u1);
    signal(2, 2) = current_signal_u2(cell_u1);
    
    power(1, 1) = signal(1, 1);
    power(1, 2) = signal(2, 2);
    
    for cell_u2  = 1:7
        if cell_u2 ~= cell_u1
            
            signal(2, 1) = current_signal_u2(cell_u2);
            signal(1, 2) = current_signal_u1(cell_u2);
            
            power(1, 2) = power(1, 2) + signal(2, 1);
            power(1, 1) = power(1, 1) + signal(1, 2);
            
            cell_u2_new = cell_u2;
        else
            cell_u2_new = 8;
        end
        
       %% Calculate alpha and floor power and alpha
        [ alpha_floor, ~ ] = calculate_alpha( signal(:,:) ); % (0 to 1 incremented by 0.1)
        power_floor = zeros(1, 2);  % (-10 to 30 incremented by 1)
        
        for i = 1:2
            
            % calculate outer cell
            power_macro = 0;
            for macro = 1:numel( all_signal_power_outer(:, 1, 1) )
                % choose randomly 2 out of 7 cells
                a = randi([1 7]);
                b = randi([1 7]);
                while a~=b
                    b = randi([1 7]);
                end
                power_macro = power_macro + all_signal_power_outer(macro, combination(i), a) + all_signal_power_outer(macro,  combination(i), b);
            end
            
            power_floor(i) = floor( 10*log10( power(i) / ( 10^( noise / 10 ) + power_macro) ) );
            if power_floor(i) >= 30
                power_floor(i) = 30;
            elseif power_floor(i) <= -10
                power_floor(i) = -10;
            end
        end
        
       %% Find the best modulation
        [modulation, throughput] = find_best_mod( power_floor, alpha_floor, 0, ccc_table );
        [modulation_jd, throughput_jd] = find_best_mod( power_floor, alpha_floor, 1, ccc_table );
        
        if throughput > best_throughput
            best_throughput = throughput;
            best_mod = modulation;
            cell_pair = [ cell_u1 cell_u2_new ];
            best_pow = power_floor;
            best_alpha = alpha_floor;
        end
        
        if throughput_jd > best_throughput_jd
            best_throughput_jd = throughput_jd;
            best_mod_jd = modulation_jd;
            cell_pair_jd = [ cell_u1 cell_u2_new ];
            best_pow_jd = power_floor;
            best_alpha_jd = alpha_floor;
        end
    end
end

% look up CCC output
ccc_output(1, 1) = ccc_table.CCCtable_conv_SINRp_alphap_QAMq_QAMp( best_pow(1) + 11, round(10*(1-best_alpha(1))) + 1, best_mod(2), best_mod(1));
ccc_output(1, 2) = ccc_table.CCCtable_conv_SINRp_alphap_QAMq_QAMp( best_pow(2) + 11, round(10*(1-best_alpha(2))) + 1, best_mod(1), best_mod(2));

ccc_output_jd(1, 1) = ccc_table.CCCtable_prop_SINRp_alphap_QAMq_QAMp( best_pow_jd(1) + 11, round(10*(1-best_alpha_jd(1))) + 1, best_mod_jd(2), best_mod_jd(1));
ccc_output_jd(1, 2) = ccc_table.CCCtable_prop_SINRp_alphap_QAMq_QAMp( best_pow_jd(2) + 11, round(10*(1-best_alpha_jd(2))) + 1, best_mod_jd(1), best_mod_jd(2));

connection(1, u1) = cell_pair(1);
connection(1, u2) = cell_pair(2);

connection_jd(1, u1) = cell_pair_jd(1);
connection_jd(1, u2) = cell_pair_jd(2);
best_mod
end

