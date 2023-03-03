function generate_plots(t,y,L,setpoints,setpoints_T)
% Generate plots 

    figure    
    adjust_plot()

    c_seaGreen = [43, 138, 87]/255;
    c_tabBlue = [30, 117, 179]/255;
    
    %Plot trajectory of photophilic strain ratio
    yyaxis left
    plot(t/60,y(:,end),'-','Color',c_seaGreen,'Linewidth',6,'MarkerSize',10)
    hold on
    
    
    %Plot the target setpoint 
    if exist('setpoints','var')
        if ~exist('setpoint_T','var') && length(setpoints) == 1
            %Single target setpoint 
            yline(setpoints,'--','Linewidth',4,'Color',[0.25, 0.25, 0.25],'FontSize',50)%,'Setpoint'
        elseif ~exist('setpoints_T','var') && length(setpoints) > 1
            error("If there are multiple setpoints, provide also the vector of timepoints when the setpoints change.")
        else
            %Dynamically changing target setpoints
            setpoint_time_intervals = {};
            for i=1:length(setpoints_T)-1
                setpoint_time_intervals{i} = [setpoints_T(i),setpoints_T(i+1)]/60;
            end
            setpoint_time_intervals{length(setpoints_T)} = [setpoints_T(length(setpoints_T)),t(end)]/60;
            for i=1:length(setpoints)
                plot(setpoint_time_intervals{i},[setpoints(i),setpoints(i)],'--','Color',[0.25, 0.25, 0.25],'Linewidth',4,'MarkerSize',10)
            end
        end
    end

    ylabel('Photophilic Fraction','Interpreter','latex')
    ylim([-0.05,1.05])
    yticks(linspace(0,1,11));
    grid on
    set(gca,'ycolor',c_seaGreen)
    set(gca,'FontSize',50)

    %Plot the delivered light inputs
    yyaxis right
    stairs(t/60,L,'-','Linewidth',4,'Color',c_tabBlue)   
    xlabel('Time [h]','Interpreter','latex')% [1/min]
    ylabel('Blue Light [a.u.]','Interpreter','latex')
    
    ylim([-40,840])
    yticks(linspace(0,800,11));
    set(gca,'ycolor',c_tabBlue) 
    set(gca,'FontSize',50)

    
    %Plot also the dynamics of resistance levels
    c_Red = [255, 42, 42]/255;
    figure
    adjust_plot()
    yyaxis left
    plot(t/60,y(:,3),'-','Color',c_Red,'Linewidth',6,'MarkerSize',10)
    ylabel('CAT [nM]','Interpreter','latex')
    ylim([0-0.05*1200,1200+0.05*1200])
    yticks(linspace(0,1200,9));
    grid on
    set(gca,'ycolor',c_Red)
    set(gca,'FontSize',50)

    yyaxis right
    stairs(t/60,L,'-','Linewidth',4,'Color',c_tabBlue)   
    xlabel('Time [h]','Interpreter','latex')% [1/min]
    ylabel('Blue Light [a.u.]','Interpreter','latex')
    
    ylim([-40,840])
    yticks(linspace(0,800,11));
    set(gca,'ycolor',c_tabBlue) 
    set(gca,'FontSize',50)
end

