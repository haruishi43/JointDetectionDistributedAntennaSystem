function modulation = find_best_mod( power, alpha, jd, ccc_table )
% Find the best modulation using CCC lookup table (ccc_table)

modulation = zeros(1, 2);

if jd == 0
    mod_list1 = squeeze(ccc_table.CCCtable_conv_SINRp_alphap_QAMq_QAMp( power(1) + 11, round( 10*(1-alpha(1)) ) + 1, :, :));
    mod_list2 = squeeze(ccc_table.CCCtable_conv_SINRp_alphap_QAMq_QAMp( power(2) + 11, round( 10*(1-alpha(2)) ) + 1, :, :));
    
    tot_mod_list = mod_list1 + mod_list2';
    [~, index] = max( tot_mod_list(:) );
    [row, col] = ind2sub( size(tot_mod_list), index );
    
    modulation(1, :) = [ col row ];
    
elseif jd == 1
    mod_list1 = squeeze(ccc_table.CCCtable_prop_SINRp_alphap_QAMq_QAMp( power(1) + 11, round( 10*(1-alpha(1)) ) + 1, :, :));
    mod_list2 = squeeze(ccc_table.CCCtable_prop_SINRp_alphap_QAMq_QAMp( power(2) + 11, round( 10*(1-alpha(2)) ) + 1, :, :));
    
    tot_mod_list = mod_list1 + mod_list2';
    [~, index] = max( tot_mod_list(:) );
    [row, col] = ind2sub( size(tot_mod_list), index );
    
    modulation(1, :) = [ col row ];
    
end

end

