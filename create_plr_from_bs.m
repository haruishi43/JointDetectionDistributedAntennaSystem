function plr_from_bs = create_plr_from_bs( antenna_coordinates, user_coordinates )
% PLR (Packet Loss Ratio)
%
% Packet loss occurs when one or more packets of data travelling across a
% network fails to reach their destination.
% 
% TODO: add equation discriptions
% 

if nargin < 2
    error('Error: not enough inputs.');
end

num_user = numel(user_coordinates);
num_cell = numel(antenna_coordinates);

distance_from_bs = zeros(num_cell, num_user);
plr_from_bs = zeros(num_cell, num_user);

% FIXME: remove magic numbers
for i = 1:num_cell
    
    distance = abs( user_coordinates - repmat( antenna_coordinates(i), 1, num_user ).' );
    
    distance_from_bs(i, :) = abs( sqrt(distance(:).^2 + 8.5^2) );
    
    plr_from_bs(i, :) = 140.7 + 36.7 * log10( distance_from_bs(i, :) * 0.001 );
    
end

end

