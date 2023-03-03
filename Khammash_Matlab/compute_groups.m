function [g_ON,phi_S,phi_R,gr,h_ON] = compute_groups(y,p_var,p_fix,L)
% This function computes useful quantities for several modeling files:
% - g_ON: fraction of CAT gene copies that are bound by active polymerase
% - phi_S: synthetic proteome fraction
% - phi_R: ribosomal proteome fraction
% - gr: growth rate
% - h_ON: light-dependent dimerization rate of the opto-T7 monomers

% If computation at a single timepoint is required, the state vector must
% be a vector.
% For computation of a time-series, provide a matrix where the columns are
% the species and the rows are the time points
    s = size(y);

    if s(2)==1
        y = y';
    end
    
    %Unpack required parameters
%     a_T = p_var(1);
%     a_C = p_var(2);
    K_G = p_var(3);
    hon_min = p_var(4);
    hon_max = p_var(5);
    K_L_on = p_var(6);
    n_L_on = p_var(7);
    n_G = p_var(8);
    h_C = p_var(9);
%     mol2fluo = p_var(10);
%     L_0 = p_var(11);
    
    N_p = p_fix(1);
    K_D = p_fix(2);
    K_C = p_fix(3);
    k = p_fix(4);
%     n_r = p_fix(5);
    n_T = p_fix(6);
    n_C = p_fix(7);
    g_0 = p_fix(8);
    phi_max = p_fix(9);
    phi_0 = p_fix(10);
    d_cell = p_fix(11);
    nu = p_fix(12);
    A_E = p_fix(13);

    % State vector
    T_T = y(:,1);
    T_D = y(:,2);
    C = y(:,3);
    
    g_ON = N_p.*T_D.^n_G./(T_D.^n_G+K_G.^n_G);
    phi_S = (n_T.*T_T + n_C.*C)./d_cell;
    gr = (phi_max - phi_0 - phi_S).*nu.*g_0./(g_0 + nu.*(1 + (A_E./K_D)./(1 + (C./(k*K_C)).^h_C)));
    phi_R = phi_0 + gr./g_0.*(1 + (A_E./K_D)./(1 + (C./(k*K_C)).^h_C));
    
    
    if exist('L','var')
        s_L = size(L);
        if s_L(2)==1
            L = L';
        end
        h_ON = hon_min + (hon_max - hon_min).*L.^n_L_on./(L.^n_L_on + K_L_on.^n_L_on);
        if isscalar(h_ON)
            h_ON = ones(length(gr),1)*h_ON;
        end
    else
        h_ON = zeros(length(gr),1);
    end
end

