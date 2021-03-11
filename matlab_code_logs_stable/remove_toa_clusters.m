function cleaned_toas = remove_toa_clusters(toas, z_thresh)
    cleaned_toas(1,:) = toas(1,:);
    toas_idx = 2;
    toas_idx_prev = 1;
    cleaned_toas_idx = 2;
    while(toas_idx < size(toas, 1))
        while(((toas(toas_idx,1) - toas(toas_idx_prev,1))*86400) < z_thresh)
            toas_idx = toas_idx + 1;
            if(toas_idx == size(toas, 1))
                return
            end
        end
        toas_idx_prev = toas_idx;
        cleaned_toas(cleaned_toas_idx, :) = toas(toas_idx, :);
        cleaned_toas_idx = cleaned_toas_idx + 1;
    end
end
      
    