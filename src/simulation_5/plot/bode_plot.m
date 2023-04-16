function bode_plot(G, legends, title, plot_name)
% This function plots the bode diagrams for magnitude and phase
% G is a vector of transfer functions
% legend is the vector of legends corresponding to the TFs G

parameters

fig = figure('Position', get(0, 'Screensize'), 'Color','w');

for i=1:length(G)
  [magG, phaseG, woutG] = bode(G(i));

  % Magnitude
  subplot(2,1,1) 
  semilogx(woutG, 20*log10(squeeze(magG)), 'LineWidth', line_width)
  hold on
  if i==length(G)
    xlabel('$\omega$ [rad/s]', 'FontSize', font_size, 'interpreter','latex')
    ylabel('Mag. [dB]', 'FontSize', font_size, 'interpreter','latex')
    grid on
    set(gca, 'FontSize', font_size)
    legend(legends,'interpreter','latex','FontSize', font_size, 'location', 'best')
  end
  
  % Phase
  subplot(2,1,2)
  semilogx(woutG, squeeze(phaseG), 'LineWidth', line_width)
  hold on
  if i == length(G)
    xlabel('$\omega$ \ [rad/s]', 'FontSize', font_size, 'interpreter','latex')
    ylabel('Phase [deg.]', 'FontSize', font_size, 'interpreter','latex')
    grid on
    set(gca, 'FontSize', font_size)
    sgtitle(title, 'FontSize', font_size);
  end

  if simulation.print_figure == 1
    fig_name = strcat(path_images, '\', plot_name,'.svg');
    export_fig(fig, fig_name);
  end

end