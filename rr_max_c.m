function [ signal, power, connection ] = rr_max_c( num_users, combination, current_signal )
% RR user scheduling with Max-C (not Max-C/I) scheduling for DA selection
% 2 users are seleted in this RR scheduling
%
% returns signal (pair of mai signal and interference), 
% cumulative power (S + I), and connections
%

% user 1 is always prioritized than user 2
u1 = combination(1);
u2 = combination(2);

connection = 8 * ones(1, num_users);
signal = zeros(2, 2); % 1 is main signal, 2 is interference
% [1, 1] is user 1's main signal, [1, 2] is user 1's interference
% [2, 1] is user 2's main signal, [2, 2] is user 2's interference
power = zeros(1, 2);

% max-c (best signal power is chosen)
current_signal_u1 = current_signal(u1, : );
[s, i] = max( current_signal_u1 );

connection(1, u1) = i;
signal(1, 1) = s;
signal(2, 2) = current_signal(u2, i);

% power (in dB)
power(1, 1) = 10*log10( abs(s) );
power(1, 2) = 10*log10( abs(signal(2, 2)) );

% max-c for user 2
current_signal_u2 = current_signal(u2, :);
[s, i] = max( current_signal_u2 );

% if user 2 tries to connect to user 1's DA, it would not let it happen
if isempty( find( i == connection(1, :), 1 ) ) == 1
    connection(1, 2) = i;
    
    signal(2, 1) = s;
    signal(1, 2) = current_signal(u1, i);
    
    power(1, 2) = power(1, 2) + 10*log10( abs(s) );
    power(1, 1) = power(1, 1) + 10*log10( abs(signal(1, 2)) );
end

end

