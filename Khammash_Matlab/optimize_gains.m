%% Parameter definition

%Load parameter values for the ODE model of the co-culture
define_parameters;

%growth rate of constitutive strain
% bJAG367:  
% gr_constitutive = 1.21424238/60;
% bJAG235:  
gr_constitutive = 1.58369475/60;
% bJAG237: 
% gr_constitutive = 1.12677044/60;

p_fix = [p_fix;gr_constitutive];


%Boolean: true if the optimal solution should be plotted after optimization
plot_opt_sol = true;

% Set point
p_cont.y_set = [0.2,0.8];

% Sampling time [min]
p_cont.dt = 30;

% Initial ratio (inoculation ratio)
f_0 = 0.5;

% Time of simulation [min]
t_final = 30*60;

% Pack everything into single parameter structure
params.p_var = p_var;
params.p_fix = p_fix;
params.p_cont = p_cont;
params.t_final = t_final;
params.f_0 = f_0;

%% Optimization
%Number of random samples from parameter space
num_particles = 50000;

%Output file
out_file = 'opt_params_cont.mat';

%Optimization
% Check if optimized variables exist in folder. If so, take the
% residual norm of the existing optimal parameters as best existing guess. 
% If not, set it to infinity.
if isfile(out_file)
    % File exists.
    load(out_file,'resnorm_best')
    resnorm_best
else
    % File does not exist.
    resnorm_best = inf
end   

% Initial guess
%Proportional gain
p_cont.k_p = 1.5327e+03;
%Integral gain
p_cont.k_i = 9.6743;
%Derivative gain
p_cont.k_d = 9.5689e+04;
%Backcalculation gain
p_cont.k_bc = 1e-1;

x0 = [p_cont.k_p,p_cont.k_i,p_cont.k_d,p_cont.k_bc];

%Parameter bounds
lb = [1e-5,...
      1e-5,...
      1e-5,...
      1e-2];  
ub = [1e7,...
      1e7,...
      1e7,...
      1e0];
  
options = optimoptions('lsqnonlin','FiniteDifferenceStepSize',100);

%First, try out the gains of the initial guess
try
    [x,resnorm_init,residual,exitflag,output,lambda,jacobian] = ...
    lsqnonlin(@(x)min_fun(x,params),x0,lb,ub,options);
catch
    resnorm_init = inf;
end

if resnorm_init < resnorm_best
    resnorm_best = resnorm_init
    save(out_file);
end

%Second, perform random sampling
tic
for i = 1:num_particles
    i
    %Update Initial guess
    random_exponents = (log10(ub)-log10(lb)).*rand(1,4) + log10(lb);
    x0 = 10.^random_exponents;
 
    %Optimization
    try
        [x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
        lsqnonlin(@(x)min_fun(x,params),x0,lb,ub,options);
        resnorm
    catch
        resnorm = inf
    end
    t_opt = toc
    if resnorm < resnorm_best
        resnorm_best = resnorm
        save(out_file);
    end
end

%% If required, simulate and plot optimal solution
if plot_opt_sol

    closed_loop = true;
    
    x = load(out_file,'x').x;
    
    y_sets = p_cont.y_set;
    
    for i=1:length(y_sets)
        %Update controller object for current target setpoint
        p_cont_opt = toggleStructArray_P_cont(x,y_sets(i),p_cont.dt);
        [t_opt,y_opt,L_opt] = simulate_timeCourse(p_var, p_fix, p_cont_opt, t_final, f_0, closed_loop,y_sets(i));

        generate_plots(t_opt,y_opt,L_opt)
    end
end
