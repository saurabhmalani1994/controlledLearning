function y0 = get_init_cond(p_var,p_fix)
% Simulate initial steady-state (assumed to be a steady-state corresponding
% to the ambient light value L0 of the pre-culture condition)

% Hypothetical light intensity of the pre-culture (cells incubated at 37 
% degrees with shaking in an incubator with transparent lid, exposed to 
% ambient light)
    L0 = p_var(end);
    
% Starting condition with all zeros     
    options = odeset('NonNegative',1);
    zero_state = zeros(1,3);
    
    t_end = 100000;
    
    tspan = linspace(0,t_end,10);
    
    lastwarn("");
    
    [~,y] = ode15s(@(t,y) ODEs_Internal(t,y,p_var,p_fix,L0), tspan, zero_state, options);
    
    if ~isempty(lastwarn)
        lastwarn("");
        [~,y] = ode23s(@(t,y) ODEs_Internal(t,y,p_var,p_fix,L0), tspan, zero_state);
    end
    
    if ~isempty(lastwarn)
        error("ode23s caused warning for initial steady-state")
    end
    
    y0 = y(end,:);
end