function coordinates = create_bs_coordinate( optional_intersite_distance )
% Function that creates the base station coordinates
%
%   No inputs are necessary when using the default value for inter-site
%   distance.
%   The number of cell is 7.
%
%   Default values:
%     Inter-site Distance: 500
%

if nargin > 0
    intersite_distance = optional_intersite_distance;
else
    intersite_distance = 500;
end

no_cell = 7; 
coordinates = zeros(no_cell, 1);

coordinates(1) = 0;         % cell 1 is the center

for a = 2:7
    coordinates(a) = intersite_distance * cos(a * pi/3 - pi/6) ...
                             + 1i * intersite_distance * sin(a * pi/3 - pi/6);
end

end

