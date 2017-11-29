function [list, tot] = create_combination_table(max,r)
% Creates a list of all the permutations given size and r

if nargin < 2 
    error('Error: not enough inputs.');
end

a = nchoosek(1:max, r);
list = zeros(0, r);

for i=1:size(a,1)
    pi = perms(a(i,:));
    list = unique([list; pi], 'rows');
end

tot = numel(list(:, 1));

end

