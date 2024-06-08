[f,xi] = ksdensity(z(z~=0));
close all
histogram(z(z~=0), 50, 'Normalization', 'pdf')
hold on
plot(xi,f, 'LineWidth',2,'color','black')
set(gca, 'Fontsize', 16)

if save_figs == 1
saveas(gcf,strcat('_results/shock_histogram'),'png')
end