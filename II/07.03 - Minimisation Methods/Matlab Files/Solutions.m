format long

syms x y z
f4(x,y) = x + y + x^2/4 - y^2 + (y^2-x/2)^2;
f5(x,y) = (1-x)^2 + 80*(y-x^2)^2;
f6(x,y,z) = 0.4*x^2 + 0.2*y^2 + z^2 + x*z;

plot_graphs = true;

% Question 1

if plot_graphs
    f = figure(1);
    f.Position = [20, 1, 1000, 500];
    f.Visible = 'off';
    X = linspace(-1.5, 1.5, 100);
    Y = linspace(-1.5, 1.5, 100);
    [X,Y] = meshgrid(X,Y);
    Z4 = double(f4(X,Y));
    subplot(1,2,1);
    contourf(X,Y,Z4,20);
    daspect([1, 1, 1])
    xlabel('x')
    ylabel('y')
    colorbar

    subplot(1,2,2);
    s = surf(X,Y,Z4);
    s.EdgeColor = 'none';
    xlim([-1.5 1.5])
    ylim([-1.5 1.5])
    xticks(-1.5:0.5:1.5)
    yticks(-1.5:0.5:1.5)
    daspect([1, 1, 3])
    xlabel('x')
    ylabel('y')
    zlabel('f(x,y)')

    print("contour1","-depsc");
    
    
    clf(f,'reset')
    f = figure(1);
    f.Position = [20, 350, 1000, 500];
    f.Visible = 'off';
    Z5 = double(f5(X,Y));
    subplot(1,2,1);
    [~, h] = contourf(X,Y,Z5,20);
    h.LevelList = [0 0.2 5 h.LevelList];
    daspect([1, 1, 1])
    xlabel('x')
    ylabel('y')
    colorbar

    subplot(1,2,2);
    s = surf(X,Y,Z5);
    s.EdgeColor = 'none';
    xlim([-1.5 1.5])
    ylim([-1.5 1.5])
    xticks(-1.5:0.5:1.5)
    yticks(-1.5:0.5:1.5)
    daspect([1, 1, 350])
    xlabel('x')
    ylabel('y')
    zlabel('f(x,y)')
    print("contour2","-depsc");
end

% Question 2

[xs1, fs1] = SD(f4, [-1.0 -1.3], 'auto', true,...
    'maxIter', 10, 'printIteration', true,...
    'functionTolerance', 0, 'stationaryTolerance', 0);

if plot_graphs
    clf('reset')
    f = figure(1);
    f.Visible = 'off';
    X = linspace(-1.1, -0.2, 100);
    Y = linspace(-1.4, -0.5, 100);
    [X,Y] = meshgrid(X,Y);
    Z4 = double(f4(X,Y));
    contourf(X,Y,Z4,20);
    hold on
    [~, h] = contour(X,Y,Z4,[fs1(end) fs1(end)]);
    h.LineWidth = 1;
    h.LineColor = [0.75, 0.75, 0.75];
    hold off
    line(xs1(:,1),xs1(:,2), 'color', 'white', 'lineWidth', 2)
    xlim([-1.1, -0.2]);
    ylim([-1.4, -0.5]);
    daspect([1, 1, 1])
    xlabel('x')
    ylabel('y')
    colorbar
    print("contour3","-depsc");
    
    
    clf('reset')
    f = figure(1);
    f.Visible = 'off';
    plot(0:1:10,fs1)
    grid on
    daspect([1 1 1])
    ylim([-1.5 1.5])
    xlim([0 10])
    yticks(-1.5:0.5:1.5)
    line([-10 10], [-1.1 -1.1], 'color', 'black')
    xlabel('n')
    ylabel('$f_4\left(x_n\right)$', 'Interpreter', 'latex')
    print("flimit1","-depsc");
end

% Question 3

[xs2, ~] = SD(f5, [0.676 0.443], 'auto', true,...
    'functionTolerance', 1e-6, 'precision', 5e-3);

if plot_graphs
    clf('reset')
    f = figure(1);
    f.Visible = 'off';
    X = linspace(0.4, 1.4, 100);
    Y = linspace(0.4, 1.4, 100);
    [X,Y] = meshgrid(X,Y);
    Z5 = double(f5(X,Y));
    [~, h] = contourf(X,Y,Z5,20);
    h.LevelList = [0 0.2 5 h.LevelList(2:end)];
    line(xs2(:,1),xs2(:,2), 'color', 'white', 'lineWidth', 2)
    daspect([1, 1, 1])
    xlabel('x')
    ylabel('y')
    colorbar
    print("contour4","-depsc");
end

% Question 4

[xs3, fs3] = CG(f4, [-1.0 -1.3], 'auto', true,...
    'maxIter', 10, 'printIteration', true,...
    'functionTolerance', 0, 'stationaryTolerance', 0);

if plot_graphs
    clf('reset')
    f = figure(1);
    f.Visible = 'off';
    X = linspace(-1.1, -0.2, 100);
    Y = linspace(-1.4, -0.5, 100);
    [X,Y] = meshgrid(X,Y);
    Z4 = double(f4(X,Y));
    contourf(X,Y,Z4,20);
    hold on
    [~, h] = contour(X,Y,Z4,[fs3(end) fs3(end)]);
    h.LineWidth = 1;
    h.LineColor = [0.75, 0.75, 0.75];
    hold off
    line(xs3(:,1),xs3(:,2), 'color', 'white', 'lineWidth', 2)
    xlim([-1.1, -0.2]);
    ylim([-1.4, -0.5]);
    daspect([1, 1, 1])
    xlabel('x')
    ylabel('y')
    colorbar
    print("contour5","-depsc");
    
    
    clf('reset')
    f = figure(1);
    f.Visible = 'off';
    plot(0:1:10,fs3)
    grid on
    daspect([1 1 1])
    ylim([-1.5 1.5])
    xlim([0 10])
    yticks(-1.5:0.5:1.5)
    line([-10 10], [-1.09527539 -1.09527539], 'color', 'black')
    xlabel('n')
    ylabel('$f_5\left(x_n\right)$', 'Interpreter', 'latex')
    print("flimit2","-depsc");
end

% Question 5

[xs4, ~] = CG(f5, [0.676 0.443], 'auto', true,...
    'functionTolerance', 1e-6, 'precision', 5e-3);

if plot_graphs
    clf('reset')
    f = figure(1);
    f.Visible = 'off';
    X = linspace(0.4, 1.4, 100);
    Y = linspace(0.4, 1.4, 100);
    [X,Y] = meshgrid(X,Y);
    Z5 = double(f5(X,Y));
    [~, h] = contourf(X,Y,Z5,20);
    h.LevelList = [0 0.2 5 h.LevelList(2:end)];
    line(xs4(:,1),xs4(:,2), 'color', 'white', 'lineWidth', 2)
    daspect([1, 1, 1])
    xlabel('x')
    ylabel('y')
    colorbar
    print("contour6","-depsc");
end

% Question 6

% [xs, fs] = DFP(f6, [1 1 1], 'auto', false, 'printIteration', true, 'printH', true, 'maxLambda', 5);

% Question 7

[xs5, fs5] = DFP(f4, [-1.0 -1.3], 'auto', true,...
    'maxIter', 10, 'printIteration', true,...
    'functionTolerance', 0, 'stationaryTolerance', 0,...
    'printH', true);

if plot_graphs
    clf('reset')
    f = figure(1);
    f.Visible = 'off';
    X = linspace(-1.1, -0.2, 100);
    Y = linspace(-1.4, -0.5, 100);
    [X,Y] = meshgrid(X,Y);
    Z4 = double(f4(X,Y));
    contourf(X,Y,Z4,20);
    hold on
    [~, h] = contour(X,Y,Z4,[fs5(end) fs5(end)]);
    h.LineWidth = 1;
    h.LineColor = [0.75, 0.75, 0.75];
    hold off
    line(xs5(:,1),xs5(:,2), 'color', 'white', 'lineWidth', 2)
    xlim([-1.1, -0.2]);
    ylim([-1.4, -0.5]);
    daspect([1, 1, 1])
    xlabel('x')
    ylabel('y')
    colorbar
    print("contour7","-depsc");
    
    
    clf('reset')
    f = figure(1);
    f.Visible = 'off';
    plot(0:1:10,fs5)
    grid on
    daspect([1 1 1])
    ylim([-1.5 1.5])
    xlim([0 10])
    yticks(-1.5:0.5:1.5)
    line([-10 10], [-1.095275394488075 -1.095275394488075], 'color', 'black')
    xlabel('n')
    ylabel('$f_6\left(x_n\right)$', 'Interpreter', 'latex')
    print("flimit3","-depsc");
end

% Question 8

[xs6, ~] = DFP(f5, [0.676 0.443], 'auto', true,...
    'functionTolerance', 1e-6, 'precision', 5e-3,...
    'printH', true);

if plot_graphs
    clf('reset')
    f = figure(1);
    f.Visible = 'off';
    X = linspace(0.4, 1.4, 100);
    Y = linspace(0.4, 1.4, 100);
    [X,Y] = meshgrid(X,Y);
    Z5 = double(f5(X,Y));
    [~, h] = contourf(X,Y,Z5,20);
    h.LevelList = [0 0.2 5 h.LevelList(2:end)];
    line(xs6(:,1),xs6(:,2), 'color', 'white', 'lineWidth', 2)
    daspect([1, 1, 1])
    xlabel('x')
    ylabel('y')
    colorbar
    print("contour8","-depsc");
end