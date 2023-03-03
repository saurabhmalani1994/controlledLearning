function adjust_plot(gray_axis)
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 28;
    grid on
    set(ax,'FontName','mwa_cmr10')
    
    if ~exist('gray_axis','var')
        gray_axis = 1;
    end
    
    if gray_axis
        set(ax, 'XColor',	[0.25, 0.25, 0.25], 'YColor',	[0.25, 0.25, 0.25])
    end
end

