function [int_val_new,L] = update_controller_directPID(int_val,y,p_cont)
% Implement PID control law to determine next light input level
    e = zeros(1,2);
    
    k_p = p_cont.k_p;
    k_i = p_cont.k_i;
    k_d = p_cont.k_d;
    k_bc = p_cont.k_bc;
    y_set = p_cont.y_set;
    dt = p_cont.dt;
    
    %Handle first time-points
    switch length(y)
        case 1
            e(1) = y_set - y(end);
            e(2) = e(1);
        otherwise
            e(1) = y_set - y(end);
            e(2) = y_set - y(end-1);
    end

    %Compute output
    int_val_new = k_i*dt*e(1) + int_val;
    ori_L = round(k_p*e(1) + int_val_new + k_d*(e(1) - e(2))/dt);   
    %Saturation
    if ori_L > 800
        L = 800;
    elseif ori_L < 0
        L = 0;
    else
        L = ori_L;
    end
    %Back-calculation (Anti-windup)
    bc = L - ori_L;
    int_val_new = int_val_new + k_bc*k_i*bc;
end

