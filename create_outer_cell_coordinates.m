function coordinates = create_outer_cell_coordinates( optional_intersite_distance )
% Create outer macro cell coordinates 

if nargin > 0
    is_d = optional_intersite_distance;
else
    is_d = 189; % inter-site distance (500m -> 189m)
end 

coordinates = zeros(6, 7);  
r = (is_d/2) / cos( pi/6 );    % cell radius

% coordinates counter-clock wise from top left
origin_points = [ -3*r + 2*is_d*1i, ...
                  -(9/2)*r - (1/2)*is_d*1i, ...
                  -(3/2)*r - (5/2)*is_d*1i, ...
                  3*r - 2*is_d*1i, ...
                  (9/2)*r + (1/2)*is_d*1i, ...
                  (3/2)*r + (5/2)*is_d*1i];

for i = 1:numel(origin_points)
    
    origin = origin_points( i );
    coordinates(i, 1) = origin;
    
    for j = 2:7
        coordinates(i, j) = origin + is_d*cos( j * pi/3 - pi/6 ) ...
                             + is_d*sin( j * pi/3 - pi/6 )*1i;
    end
    
end

end

