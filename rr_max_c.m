function [ ccc_output, ccc_output_jd, power_floor, alpha_floor, connection ] = rr_max_c( num_users, combination, current_signal, all_signal_power_outer, noise, ccc_table )
% RR user scheduling with Max-C (not Max-C/I) scheduling for DA selection
% 2 users are seleted in this RR scheduling
%
% returns power (floored) and alpha (floored)
%

%% Outputs
ccc_output = zeros(1, 2);
ccc_output_jd = zeros(1, 2);

%% RR and Max-C
% user 1 is always prioritized than user 2
u1 = combination(1);
u2 = combination(2);

connection = 8 * ones(1, num_users);
signal = zeros(2, 2); % 1 is main signal, 2 is interference
% [1, 1] is user 1's main signal, [1, 2] is user 1's interference
% [2, 1] is user 2's main signal, [2, 2] is user 2's interference
power_real = zeros(1, 2);

% max-c (best signal power is chosen)
current_signal_u1 = current_signal(u1, : );
[s, i] = max( current_signal_u1 );

connection(1, u1) = i;
signal(1, 1) = s;
signal(2, 2) = current_signal(u2, i);

% power (in dB)
power_real(1, 1) = s;
power_real(1, 2) = signal(2, 2);

% max-c for user 2
current_signal_u2 = current_signal(u2, :);
[s, i] = max( current_signal_u2 );

% if user 2 tries to connect to user 1's DA, it would not let it happen
if isempty( find( i == connection(1, :), 1 ) ) == 1
    connection(1, u2) = i;
    
    signal(2, 1) = s;
    signal(1, 2) = current_signal(u1, i);
    
    power_real(1, 2) = power_real(1, 2) + s;
    power_real(1, 1) = power_real(1, 1) + signal(1, 2);
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
    
    power_floor(i) = floor( 10*log10( power_real(i) / ( 10^( noise / 10 ) + power_macro) ) );
    if power_floor(i) >= 30
        power_floor(i) = 30;
    elseif power_floor(i) <= -10
        power_floor(i) = -10;
    end
end

%% Find the best modulation
[modulation, ~] = find_best_mod( power_floor, alpha_floor, 0, ccc_table );
[modulation_jd, ~] = find_best_mod( power_floor, alpha_floor, 1, ccc_table );

% look up CCC output
ccc_output(1, 1) = ccc_table.CCCtable_conv_SINRp_alphap_QAMq_QAMp( power_floor(1) + 11, round(10*(1-alpha_floor(1))) + 1, modulation(2), modulation(1));
ccc_output(1, 2) = ccc_table.CCCtable_conv_SINRp_alphap_QAMq_QAMp( power_floor(2) + 11, round(10*(1-alpha_floor(2))) + 1, modulation(1), modulation(2));

ccc_output_jd(1, 1) = ccc_table.CCCtable_prop_SINRp_alphap_QAMq_QAMp( power_floor(1) + 11, round(10*(1-alpha_floor(1))) + 1, modulation_jd(2), modulation_jd(1));
ccc_output_jd(1, 2) = ccc_table.CCCtable_prop_SINRp_alphap_QAMq_QAMp( power_floor(2) + 11, round(10*(1-alpha_floor(2))) + 1, modulation_jd(1), modulation_jd(2));

end

