graph_a_values([0 10], 1e-7, 1e-8, [0.5, 0.9, 1, 10], 0,... 
    'X(0) = 0', 'x0=0 apos')
graph_a_values([0 10], 1e-7, 1e-8, [-0.5, -0.9, -1, -10], 0,... 
    'X(0) = 0', 'x0=0 aneg')

graph_a_values([0 10], 1e-14, 1e-16, [0.5, 0.9, 1, 10], 0,... 
    'X(0) = 0', 'x0=0 apos high')
graph_a_values([0 10], 1e-14, 1e-16, [-0.5, -0.9, -1, -10], 0,... 
    'X(0) = 0', 'x0=0 aneg high')

graph_x_values([0 10], 1e-7, 1e-8, 0.5, linspace(0,1,5),... 
    'a=0.5', 'a=05')
graph_x_values([0 10], 1e-7, 1e-8, 0.9, linspace(0,1,5),... 
    'a=0.9', 'a=09')
graph_x_values([0 10], 1e-7, 1e-8, 1, linspace(0,1,5),... 
    'a=1', 'a=1')

table = log_graph_a_values([0 exp(12)], 1e-7, 1e-8, logspace(log10(0.01),log10(0.2),11), 0,... 
    'log(X(t)) against log(t)', 'log graphs');

hold off

x = linspace(-3,3, 300);
y1 = 2*cos(2*pi*x)-1;
y2 = 1*cos(2*pi*x)-1;
y3 = 0.5*cos(2*pi*x)-1;

figure('visible','off')
plot(x,y1)
xticks(-10:1:10)
xlabel("\chi")
yticks(-10:1:10)
yl = ylabel("$\frac{d\chi}{dt}$",'Interpreter','latex');
set(yl,'rotation',0,'VerticalAlignment','bottom')
axis([-3, 3, -4, 4])
line([-10, 10], [0, 0], 'color', 'black')
line([-10, 10], [-1, -1], 'color', [0.7, 0.7, 0.7])
line([1/6, 1/6], [10, -10], 'color', [0.25, 0.25, 0.25], 'LineStyle', '--')
line([-1/6, -1/6], [10, 0], 'color', [0.5, 0.5, 0.5], 'LineStyle', '--')
line([5/6, 5/6], [0, -10], 'color', [0.5, 0.5, 0.5], 'LineStyle', '--')
grid on
daspect([1 2 1])
annotation('arrow',[0.495 0.54],[0.65 0.65])
annotation('arrow',[0.625 0.54],[0.2 0.2])
print("large a.eps", '-depsc');

figure('visible','off')
plot(x,y2)
xticks(-10:1:10)
xlabel("\chi")
yticks(-10:1:10)
yl = ylabel("$\frac{d\chi}{dt}$",'Interpreter','latex');
set(yl,'rotation',0,'VerticalAlignment','bottom')
axis([-3, 3, -4, 4])
line([-10, 10], [0, 0], 'color', 'black')
line([-10, 10], [-1, -1], 'color', [0.7, 0.7, 0.7])
line([0, 0], [10, -10], 'color', [0.25, 0.25, 0.25], 'LineStyle', '--')
line([1, 1], [0, -10], 'color', [0.5, 0.5, 0.5], 'LineStyle', '--')
grid on
daspect([1 2 1])
annotation('arrow',[0.647 0.5175],[0.3 0.3])
print("mid a.eps", '-depsc');

figure('visible','off')
plot(x,y3)
xticks(-10:1:10)
xlabel("\chi")
yticks(-10:1:10)
yl = ylabel("$\frac{d\chi}{dt}$",'Interpreter','latex');
set(yl,'rotation',0,'VerticalAlignment','bottom')
axis([-3, 3, -4, 4])
line([-10, 10], [0, 0], 'color', 'black')
line([-10, 10], [-1, -1], 'color', [0.7, 0.7, 0.7])
grid on
daspect([1 2 1])
annotation('arrow',[0.7765 0.2585],[0.3 0.3])
print("small a.eps", '-depsc');


function result = flowODE(t, x, a)
% Contains non-dimensionalised ODE
    result = a*cos(2*pi*(x-t));
end

function result = binary_search(arr, T)
    n = size(arr);
    L = 1;
    R = n(1);
    result = floor((L+R)/2);
    while L <= R
        m = floor((L+R)/2);
        if arr(m) >= T && arr(m-1) < T
            result = m;
            break
        elseif arr(m) < T
            L = m+1;
        else
            R = m-1;
        end
    end
end

function graph_a_values(tspan, reltol, abstol, a_values, x0,...
    graph_title, filename)
% graphs solution to the ODE

    figure('visible','off')
    opts = odeset('RelTol',reltol,'AbsTol',abstol);
    miny = inf;
    maxy = -inf;
    for a = a_values
        [t, x] = ode45(@(t,x) flowODE(t,x,a),...
                       tspan,x0,opts);
        plot(t,x, 'DisplayName',join(['a = ', num2str(a)]));
        grid on;
        miny = min([min(min(x)) miny]);
        maxy = max([max(max(x)) maxy]);
        yticks(floor(miny):1:ceil(maxy));
        xticks(floor(tspan(1)):1:ceil(tspan(2)));
        daspect([1 1 1]);
        hold on;
        legend('show');
    end
    xlabel(join([num2str(tspan(1)), "t", num2str(tspan(2))],...
                " \leq "));
    ylabel('X(t)');
    legend('Location','northwest','AutoUpdate','off');
    line([0 0], ylim, 'Color','black');
    line(xlim, [0 0], 'Color','black');
    hold off;
    title(graph_title);
    print(join(filename,'.eps'), '-depsc');
end

function graph_x_values(tspan, reltol, abstol, a, x0_values,...
    graph_title, filename)
% graphs solution to the ODE

    figure('visible','off')
    opts = odeset('RelTol',reltol,'AbsTol',abstol);
    miny = inf;
    maxy = -inf;
    for x0 = x0_values
        [t, x] = ode45(@(t,x) flowODE(t,x,a),...
                       tspan,x0,opts);
        plot(t,x, 'DisplayName',join(['X(0) = ', num2str(x0)]));
        grid on;
        miny = min([min(min(x)) miny]);
        maxy = max([max(max(x)) maxy]);
        yticks(floor(miny):1:ceil(maxy));
        xticks(floor(tspan(1)):1:ceil(tspan(2)));
        daspect([1 1 1]);
        hold on;
        legend('show');
    end
    xlabel(join([num2str(tspan(1)), "t", num2str(tspan(2))],...
                " \leq "));
    ylabel('X(t)');
    legend('Location','eastoutside','AutoUpdate','off');
    line([0 0], ylim, 'Color','black');
    line(xlim, [0 0], 'Color','black');
    hold off;
    title(graph_title);
    print(join(filename,'.eps'), '-depsc');
end

function T = log_graph_a_values(tspan, reltol, abstol, a_values, x0,...
    graph_title, filename)
% graphs solution to the ODE

    figure('visible','off')
    opts = odeset('RelTol',reltol,'AbsTol',abstol);
    miny = inf;
    maxy = -inf;
    yint = [];
    for a = a_values
        [t, x] = ode45(@(t,x) flowODE(t,x,a),...
                       tspan,x0,opts);
        x = arrayfun(@(y) log(y), x);
        t = arrayfun(@(y) log(y), t);
        ind = binary_search(t, 7);
        t1 = t(ind);
        x1 = x(ind);
        t2 = t(end);
        x2 = x(end);
        c = (x1*t2-x2*t1)/(t2-t1);
        yint = [yint; [a, c, exp(c), 0.5*a^2,...
            abs(exp(c)-0.5*a^2), abs(exp(c)-0.5*a^2)/(0.5*a^2)]];
        plot(t,x, 'DisplayName',join(['a = ', num2str(a)]));
        miny = min([min(min(x)), miny]);
        maxy = max([max(max(x)), maxy]);
        daspect([1, 1, 1]);
        hold on;
        legend('show');
    end
    T = array2table(yint,'VariableNames',{'a','y_int','calc_drift_vel',...
        'half_a_squared', 'abs_error', 'rel_error'});
    disp(T)
    xlabel(join(['0', "log(t)", num2str(log(tspan(2)))],...
                " \leq "));
    ylabel('log(X(t))');
    axis([-2, log(tspan(2)), -6, 10])
    legend('Location','northwest','AutoUpdate','off');
    line([0 0], ylim, 'Color','black');
    line(xlim, [0 0], 'Color','black');
    yticks(-100:1:20);
    xticks(-100:1:ceil(log(tspan(2))));
    grid on;
    hold off;
    title(graph_title);
    print(join(filename,'.eps'), '-depsc');
end