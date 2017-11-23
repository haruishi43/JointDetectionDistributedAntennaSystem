function impossible_list = create_impossible_pairs( num_users, num_cells )
% Function that creates all the impossible(*) antenna pairs.
%
% * impossible meaning, pairs that aren't allowed
%
% the last connection (num_cell + 1) means the user doesn't connect to any
% base station.
%


%% Variables:
pair_per_user = num_cells + 1;      % total number of choices
total_combinations = pair_per_user^num_users;       % total number of combinations

%% Simulation:
impossible_list = zeros(1, total_combinations);
count = 1;
user_pairing = zeros(1,num_users);

for pair = 1:total_combinations
    
    pair_shift = pair - 1;
    
    for i = 1:num_users
        user_pairing(i) = fix( pair_shift / pair_per_user^(num_users - i) ) + 1;
        pair_shift = rem( pair_shift, pair_per_user^(num_users - i) );
    end
    
    if user_pairing(1) ~= pair_per_user
        if user_pairing(1) == user_pairing(2)
            impossible_list(count) = pair;
            count = count + 1;
            continue;
        elseif user_pairing(1) == user_pairing(3)
            impossible_list(count) = pair;
            count = count + 1;
            continue;
        end
        
    end
    
    if user_pairing(2) ~= pair_per_user
        if user_pairing(2) == user_pairing(3)
            impossible_list(count) = pair;
            count = count + 1;
            continue;
        end
    end
    
    % When every user doesn't connect to anything:
    if user_pairing(1) == pair_per_user && user_pairing(2) == pair_per_user && user_pairing(3) == pair_per_user
        impossible_list(count) = pair;
        count = count + 1;
        continue;
    end 
end

% remove unwanted zeros:
impossible_list( impossible_list == 0) = [];

end

