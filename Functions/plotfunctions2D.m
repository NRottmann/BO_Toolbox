% Creates 2D-Plots of all function within this folder and stores them in
% the folder plots

map = [  0,  75,  90; % dark blue
       149, 188,  14; % green
       250, 187,   0; % yellow
       236, 116,   4; % orange
       ]/255;
map = [map; ones(5, 3)];
map = [interp(map(:, 1), 100) interp(map(:, 2), 100) interp(map(:, 3), 100)];
map = min(max(map(1:300, :), 0), 1);

listing = dir();
for idx=1:length(listing)
    name = listing(idx).name;
    if listing(idx).isdir || startsWith(name, 'test_fun') || startsWith(name, 'plot')
        continue
    end
    
    name = name(1:end-2); % remove .m extension
    if startsWith(name, 'AirHockey')
        continue
    end
    func = str2func(name);
    if nargin(func) == 1
        f = func(2);
    else
        f= func();
    end
    vars = f.vars;
    x_range = vars(1).Range;
    y_range = vars(2).Range;
    [X,Y] = meshgrid(linspace(x_range(1), x_range(2)),...
        linspace(y_range(1), y_range(2)));
    x = X(:);
    y = Y(:);
    z = zeros(size(x));
    for i=1:length(x)
        z(i) = f.call(x(i), y(i));
    end
    Z = reshape(z, size(X));
    width = 7.4;        % centimeter for a half page site
    height = 0.75*width;
    figure('Units','centimeters','Position', [2 2 width height])
    set(gca ,'FontName', 'Latin Modern Mono Light', 'FontSize' ,10);
    surfc(X,Y,Z, 'EdgeAlpha', 0.2)
    %title(name);
    xlabel('x');
    ylabel('y');
    zlabel(strcat(lower(name), {' '}, '(x,y)'));
%     xlabel('$$x$$', 'interpreter','latex');
%     ylabel('$$y$$', 'interpreter','latex');
%     zlabel(strcat(lower(name), {' '}, '$$(x,y)$$'),'interpreter','latex')%, 'FontName', 'Latin Modern Mono Light', 'FontSize' ,10);
    colormap(map)
    saveas(gcf, fullfile('plots', name), 'png');
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf, fullfile('plots', name), '-dpdf', '-r500');
    close all
end
close all