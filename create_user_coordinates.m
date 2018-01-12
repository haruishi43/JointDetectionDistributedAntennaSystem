function coordinates = create_user_coordinates( antenna_coordinates, optional_num_users, optional_intersite_distance, optional_preset )
% Function that creates user coordinates
%
%   Input: 
%       Antenna Coordinates (must be 7 cells)
%
%   Optional:
%       1. Presets for User Placement  (default: placed randomly)
%       2. Number of Users             (default: 3)
%       3. Inter-site Distance         (default: 500m -> 189m)
%   
%   Outputs:
%       Coordinates for Each User
%

%% Randomize every simulation:
rng('shuffle');

%% Check input parameters:
if nargin == 0
    error('Error: must have atleast 1 input.');
end

if numel(antenna_coordinates) ~= 7
    error('Error: there should be 7 antenna coordinates.');
end

switch nargin
    case 1
        % default case:
        num_users = 3;
        is_distance = 189;
        user_placement = randi(7, num_users, 1);
    case 2
        num_users = optional_num_users;
        is_distance = 189;
        user_placement = randi(7, num_users, 1);
    case 3
        num_users = optional_num_users;
        is_distance = optional_intersite_distance;
        user_placement = randi(7, num_users, 1);
    case 4
        num_users = optional_num_users;
        is_distance = optional_intersite_distance;
        user_placement = optional_preset;
        
        if num_users ~= numel(user_placement)
            error("ERROR: the number of users and the number of elements in the preset doesn't match.");
        end
end

%% Simulation:

MIN_DISTANCE = 5;      % minimum distance
coordinates = zeros(num_users, 1);

for i = 1:num_users
    
    % initialize
    coordinates(i) = 0;
    
    while coordinates(i) == 0
        dx = (rand-0.5) * 2 * is_distance / sqrt(3);
        dy = (rand-0.5) * 2 * is_distance / sqrt(3);

        % 1. the cell is a hexagon with the radius of 289m 
        % 2. user is placed no closer than 10m from the base station
        if (abs(dx) - is_distance / 2 / sqrt(3)) * sqrt(3) > is_distance / 2 - abs(dy) || abs(dy) > is_distance / 2 ... 
            || abs(dx+dy*1i) < MIN_DISTANCE

            % redo the loop
            coordinates(i) = 0;
        else
            % place the user
            coordinates(i) = antenna_coordinates(user_placement(i)) + dx + dy*1i;
        end
    end
end

end

