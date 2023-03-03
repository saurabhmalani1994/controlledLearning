% Main

%% Parameter definition

%Define model parameters
define_parameters;

%growth rate of constitutive strain
% bJAG367:  
% gr_constitutive = 1.21424238/60;
% bJAG235:  
gr_constitutive = 1.58369475/60;
% bJAG237: 
% gr_constitutive = 1.12677044/60;

p_fix = [p_fix;gr_constitutive];


%Define controller parameters

%Manually define gains
% Parameters for figure with pathological example of PID control 1
% Proportional gain
% p_cont.k_p = 1e0;
% %Integral gain
% p_cont.k_i = 1e0;
% %Derivative gain
% p_cont.k_d = 1;
% %Backcalculation gain
% p_cont.k_bc = 1;
% 
% % Parameters for figure with pathological example of PID control 2
% % Proportional gain
% p_cont.k_p = 1e5;
% %Integral gain
% p_cont.k_i = 1;
% %Derivative gain
% p_cont.k_d = 1;
% %Backcalculation gain
% p_cont.k_bc = 1;

%Load gains of a previously optimized controller
% bJAG235
a = load('opt_params_cont_bJAG235.mat');

p_cont = toggleStructArray_P_cont(a.x,a.p_cont.y_set,a.p_cont.dt);

% Set point
p_cont.y_set = 0.7;

% Sampling time [min]
p_cont.dt = 30;

% Initial ratio (inoculation ratio)
f_0 = 0.47;

% Time of simulation [min]
t_final = 30*60;

%% Simulation Open Loop

L_OL = 0;
closed_loop = false;

[t_OL,y_OL,L_OL] = simulate_timeCourse(p_var, p_fix, p_cont, t_final, f_0, closed_loop, L_OL);
generate_plots(t_OL,y_OL,L_OL)

%% Simulation Closed Loop (single setpoint)

closed_loop = true;
[t_CL,y_CL,L_CL] = simulate_timeCourse(p_var, p_fix, p_cont, t_final, f_0, closed_loop);

generate_plots(t_CL,y_CL,L_CL,p_cont.y_set)

%% Simulation Closed Loop (changing setpoints)

closed_loop = true;
% setpoint_T = 60*[0,10,30];
% setpoint_Y = [0.2,0.8,0.4];
setpoint_T = 60*[0,12.5];
setpoint_Y = [0.7,0.3];
[t_CL,y_CL,L_CL] = simulate_timeCourse_ChangingSetpoints(p_var, p_fix, p_cont, t_final, f_0, setpoint_T,setpoint_Y);

generate_plots(t_CL,y_CL,L_CL,setpoint_Y,setpoint_T)