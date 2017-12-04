function channel_response_freq = add_rayleigh_fading( num_user, num_cell, opt_num_rb, opt_num_sc_in_rb, opt_band )
% Rayleigh Fading based on Jake's Model:
% FIXME: add desciptions

%% Randomize:
rng('shuffle');

%% Check Inputs:
if nargin < 3
    num_rb = 24;
    num_sc_in_rb = 12;
    band_per_rb = 180*10^3;
    band = band_per_rb * num_rb;
else
    num_rb = opt_num_rb;
    num_sc_in_rb = opt_num_sc_in_rb;
    band = opt_band;
end

%% Variables:
num_sc = num_rb * num_sc_in_rb;     % # of total subcarriers
num_paths = 6;                      % # of paths based on Jake's model
rms_delay_spread = 1.0 * 10^(-6);   % root mean square of delay spread
interval = 1 / band;                % sampling interval for Jake's model

% FIXME: using?
%pow_amp = [10^(-0.30) 10^(-0.00) 10^(-0.20) 10^(-0.6) 10^(-0.8) 10^(-1.0)];

%% Calculate Delay:
delay_profile = zeros(num_paths, 1);
for i = 1:num_paths
    delay_profile(i) = exp( - (i - 1) / ( rms_delay_spread / interval ) );
end
delay_profile = delay_profile / sum(delay_profile);

%% Add Rayleigh Fading:
channel_response_time = zeros(num_user, num_cell, num_sc);
channel_response_freq = zeros(num_user, num_cell, num_sc);
for user = 1:num_user
    for cell = 1:num_cell
        % 1 ~ # of paths:
        channel_response_time(user, cell, 1:num_paths) = (1/sqrt(2).*( randn(num_paths, 1) + 1i*randn(num_paths, 1) )).*sqrt(delay_profile);
        % # of paths +1 ~ # of subcarriers:
        channel_response_time(user, cell, (num_paths+1):num_sc) = zeros( 1, 1, (num_sc - num_paths) );
        
        % calculate frequency:
        channel_response_freq(user, cell, :) = fft( channel_response_time(user, cell, :) );
    end
end
end