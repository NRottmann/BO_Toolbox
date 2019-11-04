% Creates 2D-Plots of all function within this folder and stores them in
% the folder plots
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
    surf(X,Y,Z)
    %title(name);
    xlabel('$$x$$', 'interpreter','latex');
    ylabel('$$y$$', 'interpreter','latex');
    zlabel(strcat(lower(name), {' '}, '$$(x,y)$$'),'interpreter','latex', 'FontName', 'Latin Modern Mono Light', 'FontSize' ,10);
    saveas(gcf, fullfile('plots', name), 'png');
end
close all