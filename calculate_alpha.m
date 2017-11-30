function [ alpha_floor, alpha ] = calculate_alpha( signal )
% 

% calculate alpha
alpha = zeros(1, 2);
alpha_floor = zeros(1, 2);

for i = 1:2
    alpha(1, i) = signal(i, 1) / ( signal(i, 1) + signal(i, 2) );
    if isnan( alpha(1, i) )
        alpha(1, i) = 1;
    end
    
    alpha_floor(1, i) = round( alpha(1, i), 1 );
    if alpha_floor(1, i) > 1.0
        alpha_floor(1, i) = 1;
    end
end

end

