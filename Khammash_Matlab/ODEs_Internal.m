function dydt = ODEs_Internal(t,y,p_var,p_fix,L)
% This function returns the right-hand side of the ODEs that describe the
% growth-control circuit operating in the photophilic strain.

    %Unpack required parameters of the system
    a_T = p_var(1);
    a_C = p_var(2);
%     K_G = p_var(3);
%     hon_min = p_var(4);
%     hon_max = p_var(5);
%     K_L_on = p_var(6);
%     n_L_on = p_var(7);
%     n_G = p_var(8);
%     h_C = p_var(9);
%     mol2fluo = p_var(10);
%     L_0 = p_var(11);
    
%     N_p = p_fix(1);
%     K_D = p_fix(2);
%     K_C = p_fix(3);
%     k = p_fix(4);
    n_r = p_fix(5);
%     n_T = p_fix(6);
%     n_C = p_fix(7);
    g_0 = p_fix(8);
%     phi_max = p_fix(9);
%     phi_0 = p_fix(10);
    d_cell = p_fix(11);
    nu = p_fix(12);
%     A_ext = p_fix(13);
    
    %Unpack states
    T_T = y(1);
    T_D = y(2);
    C = y(3);
    
    %Compute useful quantities
    [g_ON,~,~,gr,h_ON] = compute_groups(y,p_var,p_fix,L);
    
    %Unbound ribosomes
    r_u = d_cell*gr/(n_r*g_0);
    
    %ODEs
    dT_T = a_T*r_u*gr/nu - gr*T_T;

    dT_D = h_ON*(T_T - 2*(T_D + g_ON))^2 - gr*(T_D + g_ON);     
    
    dC = a_C*r_u*g_ON/nu - gr*C;    
    
    dydt = [dT_T;...
            dT_D;...
            dC];
end

