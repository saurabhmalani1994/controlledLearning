function [t,y,L] = simulate_timeCourse(p_var, p_fix, p_cont, t_final, f_0, closed_loop, L_OL)
%Simulates a co-culture time-course experiment with a given initial
%photophilic-strain-fraction f_0 at time 0 until t_final. 

% If the boolean variable closed_loop is true, then the controller function 
% is called periodically every p_cont.dt minutes (sampling time) to update 
% the light input. If closed_loop is false, a value of L_OL has to be
% provided to simulate an open-loop experiment
    
%% Preliminary checks
    if ~exist('closed_loop','var')
        closed_loop = false;
    end
    
    if ~closed_loop && ~exist('L_OL','var') 
        error('if open-loop experiment is to be simulated, a value for the fixed light intensity must be provided as last argument.')
    end
    
    if closed_loop
        cont_fun = @(int_val,y)update_controller_directPID(int_val,y(:,4),p_cont);
    else
        cont_fun = @(l,y)update_controller_OL(L_OL);
    end
    
    %if t_final is not a multiple of p_cont.dt, round it to the closest
    %multple
    t_final=p_cont.dt*round(t_final/p_cont.dt);
    
%% Simulate
    options = odeset('NonNegative',1);
   
    % Simulate initial steady-state (assumed to be a steady-state corresponding
    % to the ambient light value L0 of the pre-culture condition)
    y0 = get_init_cond(p_var,p_fix);

    % Add initial strain ratio to inital conditions
    y0 = [y0, f_0];
    
    % output times of simulation: return num_sim_out simulation points 
    % for every sampling interval.
    num_sim_out = 3;
    output_interval = linspace(0,p_cont.dt,num_sim_out);
    
    t = linspace(0,t_final,t_final/(p_cont.dt/(num_sim_out-1))+1);
    
    % Simulation output
    y = zeros(length(t),length(y0));
    
    % System state at sample times
    y_sample = zeros(t_final/p_cont.dt+1,length(y0));
    
    % Controller output
    L = zeros(length(t),1);
    
    curr_int = 0;    
    i_sim = 1;
    i_sample = 1;    
    
    while i_sim < length(t)
        % Simulate sampling
        y_sample(i_sample,:) = y0;
        % Update controller output based on current strain ratio
        [curr_int, curr_L] = cont_fun(curr_int,y_sample(1:i_sample,:));
        % Simulate until next sampling time
        [~,y(i_sim:i_sim+num_sim_out-1,:)] = ode15s(@(t,y) ODEs(t,y,p_var,p_fix,curr_L), output_interval, y0, options);
        
        if ~isempty(lastwarn)
            lastwarn("");
            [~,y(i_sim:i_sim+num_sim_out-1,:)] = ode23s(@(t,y) ODEs(t,y,p_var,p_fix,curr_L), output_interval, y0);
        end
        if ~isempty(lastwarn)
            error("ode23s caused warning during integration.")
        end
        
        % Save simulation results
        y0 = y(i_sim+num_sim_out-1,:);
        L(i_sim:i_sim+num_sim_out-1) = curr_L;
        
        % Update time
        i_sim = i_sim+num_sim_out-1;        
        i_sample = i_sample + 1;
    end
end

