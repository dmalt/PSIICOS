function matrix = get_con_sim_matrix(CS, isUpper)

    matrix = zeros(10,10);
    l = 1;
    for i = 1:10
        for j = i+1:10
            matrix(i,j) = CS(l);    
            if ~isUpper
                matrix(j,i) = CS(l); 
            end

            l = l + 1;
        end
    end
