function [ output_args ] = calculate_sinr( selected_cell )
% Function that calculates SINR
% 
% Number of EIRP Base is 0.
%

%% Check:
if nargin < 1
   error('Error: input is needed'); 
end


%% Variables:
rnd = -174; % receiver noise density

num_eirp = 1;
eirp_base = zeros(1, num_eirp);

noise_pow = 1;


%% Simulation:


end

