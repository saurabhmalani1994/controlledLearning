function p_cont = toggleStructArray_P_cont(x,y_set,dt)
%If x passed as an array, packs controller parameters into structure
    try
        p_cont.k_p = x(1);
        p_cont.k_i = x(2);
        p_cont.k_d = x(3);
        p_cont.k_bc = x(4);
        p_cont.y_set = y_set;
        p_cont.dt = dt;
        
    catch
        %If x passed as a struct, packs controller parameters into array
        try
            p_cont(1) = x.k_p;
            p_cont(2) = x.k_i;
            p_cont(3) = x.k_d;
            p_cont(4) = x.k_bc;
        catch
            error('Controller parameters must be passed as either array or structure')
        end
    end            
end

