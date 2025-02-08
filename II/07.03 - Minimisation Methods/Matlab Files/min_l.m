function l = min_l(l_list, f_list)
    n = size(f_list);
    n = n(2);
    if f_list(1) <= f_list(2)
        l = l_list(2)/2; % assume l_list(1)=0, and don't want l=0 else
        % go nowhere but l_list(2) too large, so try half way point
    elseif f_list(n) <= f_list(n-1)
        l = l_list(n);
    else
        for j = 2:1:n-1
            if (f_list(j) <= f_list(j-1))&&(f_list(j) <= f_list(j+1))
                l = l_list(j);
                break
            end
        end
    end
end