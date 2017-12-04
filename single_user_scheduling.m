function ccc_output = single_user_scheduling( current_user, all_signal_power, all_signal_power_outer, noise, ccc_table )
% Scheduling with single user (Max-C)

current_signal = all_signal_power(current_user, :);


[signal, ~] = max( current_signal );

power_macro = 0;

for macro = 1:numel( all_signal_power_outer(:, 1, 1) )
    
    % choose randomly 1 out of 7 cells
    power_macro = power_macro + all_signal_power_outer(macro, current_user, randi([1 7]));
end

power = 10*log10( abs(signal) / ( 10^(noise/10) + abs(power_macro) ) );
power_floor = floor( power );
if power_floor >= 30
    power_floor = 30;
elseif power_floor <= -10
    power_floor = -10;
end

alpha_floor = 1;

mod_list = squeeze(ccc_table.CCCtable_conv_SINRp_alphap_QAMq_QAMp( power_floor + 11, round( 10*(1-alpha_floor) ) + 1, :, :));

[ccc_output, ~] = max( mod_list(:) );

end

