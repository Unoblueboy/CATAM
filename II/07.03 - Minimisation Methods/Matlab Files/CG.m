function [x_list, f_list] = CG(func, x0, varargin) % max_iter, auto
% Use the Conjugate Gradients Algorithm

    p = inputParser;
    addRequired(p,'func',@(f) isa(f, 'symfun'));
    addRequired(p,'x0',@isnumeric);
    addParameter(p,'maxIter',100);
    addParameter(p,'auto',false);
    addParameter(p,'stationaryTolerance', 1e-4, @(x) isnumeric(x) && (x>=0));
    addParameter(p,'functionTolerance', 1e-4, @(x) isnumeric(x) && (x>=0));
    addParameter(p,'precision', 1e-2, @(x) isnumeric(x) && (x>0));
    addParameter(p,'maxLambda', 2, @(x) isnumeric(x) && (x>0));
    addParameter(p,'printIteration', false, @islogical);
    parse(p, func, x0, varargin{:})
    
    auto = p.Results.auto;
    func = p.Results.func;
    func_tol = p.Results.functionTolerance;
    max_iter = p.Results.maxIter;
    stat_tol = p.Results.stationaryTolerance;
    x0 = p.Results.x0;
    prec = p.Results.precision;
    print_iter = p.Results.printIteration;
    maxl = p.Results.maxLambda;

    inputs = argnames(func);
    n = size(inputs);
    n = n(2);
    grad_sym = symfun.empty(n,0);
    grad_sym = grad_sym(1);
    for j = 1:n
        grad_sym(j) = diff(func, inputs(j));
    end
    g = symfun(grad_sym, inputs)';
    
    x_list = zeros(max_iter, n);
    x_list(1,:) = x0;
    temp_x = num2cell(x0);
    f_list = zeros(1, max_iter);
    f_list(1) = double(func(temp_x{:}));
    
    if print_iter
        disp("Iteration 0")
        disp(join(["f(x_0) = " num2str(f_list(1),'%.16g')]))
    end
    
    g_list = zeros(n,n-1);
    s_list = zeros(n,n-1);
    iteration = 1;
    while iteration <= max_iter
        last_x = x_list(iteration, :);
        temp_x = num2cell(last_x);
        
        % calculate search function
        grad = double(g(temp_x{:}));
        s = -grad;
        cycle = mod(iteration, n);
        if cycle > 1 % cycle ~= 1 or 0
            g_prev = g_list(:,cycle-1);
            s_prev = s_list(:,cycle-1);
            b = (grad'*grad)/(g_prev'*g_prev);
            s = -grad + b*s_prev;
        elseif cycle == 0
            g_prev = g_list(:,n-1);
            s_prev = s_list(:,n-1);
            b = (grad'*grad)/(g_prev'*g_prev);
            s = -grad + b*s_prev;
        end
        if cycle ~= 0
            g_list(:,cycle) = grad;
            s_list(:,cycle) = s;
        end
        
        if auto && (norm(s) < stat_tol)
            break
        end
        if print_iter
            disp(' ');
            disp(join(['Iteration ',num2str(iteration,'%.16g')]));
        end
        % get lambda l
        l_list = 0:prec:maxl;
        ordinates = cell(n,1);
        for j = 1:1:n
            ordinates{j} = last_x(j) + l_list*s(j);
        end
        f_l_list = func(ordinates{:});
        if ~auto
            % manual - display graph of l
            clf('reset')
            figure(1)
            plot(l_list, f_l_list)
            xlabel('\lambda');
            ylabel('f(x+\lambda s)')
            grid on
            % manual - user input value of l
            l = input('Enter value for lambda (inf to exit): ');
            if l == inf
                break
            end
        else
            % auto - check exit conditions
            l = min_l(l_list, f_l_list);
            if print_iter
                disp(join(['lambda is ', num2str(l,'%.16g')]))
            end
        end
        
        new_x = last_x + l*s';
        temp_x = num2cell(new_x);
        new_f = double(func(temp_x{:}));
        last_f = f_list(iteration);
        if print_iter
            disp(join(['f(x_' num2str(iteration,'%.16g')...
                ') = ' num2str(new_f,'%.16g')]))
            disp(join(['f(x_' num2str(iteration,'%.16g') ') - f(x_'...
                num2str(iteration-1,'%.16g') ') = '...
                num2str(new_f-last_f,'%.16g')]))
        end
        
        x_list(iteration+1,:) = new_x;
        f_list(iteration+1) = new_f;
        if auto
            last_f = f_list(iteration);
            if abs(new_f-last_f) < func_tol
                break
            end
        end
        
        iteration = iteration + 1;
    end
    x_list = x_list(1:iteration,:);
    f_list = f_list(1:iteration);
end