function plr = create_plr_from_outer_cell( antenna_coordinates, user_coordinates )
% Propagation loss ratio from the outer cell regions

if nargin < 2
    error('Error: not enough inputs.');
end

num_user = numel(user_coordinates);
num_macro_cell = numel(antenna_coordinates(:, 1));
num_antenna = numel(antenna_coordinates(1, :));

plr = zeros(num_macro_cell, num_user, num_antenna);

for i = 1:num_macro_cell
    
    distance_from_bs = zeros(num_antenna, num_user);
    
    for j = 1:num_antenna
        
        distance = abs( user_coordinates - repmat( antenna_coordinates(i, j), 1, num_user ).' );
        
        distance_from_bs(j, :) = abs( sqrt(distance(:).^2 + 8.5^2) );
        
        plr(i, :, j) = 140.7 + 36.7 + log10( distance_from_bs(j, :) * 0.001 );
        
    end

end

