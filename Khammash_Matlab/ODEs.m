function dydt = ODEs(t,y,p_var,p_fix,L)
%This function returns the right-hand side of the ODE system that simulates 
%co-culture dynamics.

%All equations apart from the last one describe the internal state of the
%photophilic strain. These are encapsulated within the function ODE_Internal.
%The last equation corresponds to the right-hand side of the 1D-ODE
%that governs the time evolution of the photophilic strain fraction in the 
%photophilic-constitutive co-culture. 

%The current growth rate of the photophilic strain is computed through the 
%function compute_groups.

    dydt = zeros(size(y));

    %Growth rate of constitutive strain
    gr_constitutive = p_fix(14);
    
    %Right-hand side corresponding to the internal state of the photophilic
    %strain
    dydt(1:end-1) = ODEs_Internal(t,y,p_var,p_fix,L);
    
    %Current growth rate of the photophilic strain     
    [~,~,~,gr_photophilic,~] = compute_groups(y,p_var,p_fix,L);
    
    %Right-hand side of the equation that governs the co-culture dynamics     
    dydt(end) = (1-y(end))*(gr_photophilic-gr_constitutive)*y(end);
end

