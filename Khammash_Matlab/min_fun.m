function out = min_fun(x,params)
%Minimization function for optimization of controller gains

    %Unpack parameters
    p_var = params.p_var;
    p_fix = params.p_fix;
    p_cont = params.p_cont;
    t_final = params.t_final;
    f_0 = params.f_0;
    y_sets = p_cont.y_set;
    
    %Output errors
    out = [];
    
    %For all required setpoints, simulate and compute error (the sum of
    %squares of these errors is the minimization objective)
    for i=1:length(y_sets)
        %Update controller object for current target setpoint
        p_cont_new = toggleStructArray_P_cont(x,y_sets(i),p_cont.dt);
    
        %Simulate
        closed_loop = true;
        [~,y,~] = simulate_timeCourse(p_var, p_fix, p_cont_new, t_final, f_0, closed_loop);

        %Compute errors
        out = [out; y(:,end)-p_cont_new.y_set];
    end
end

