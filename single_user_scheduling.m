function ccc_output = single_user_scheduling( current_user, all_signal_power, ccc_table )
% Scheduling with single user (Max-C)

current_signal = all_signal_power(current_user, :);

[signal, ~] = max( current_signal );

power = 10*log10( abs(signal) );
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

