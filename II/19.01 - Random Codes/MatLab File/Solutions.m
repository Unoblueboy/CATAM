format long

Question_1 = false;
Question_2 = false;
Question_3 = true;
Question_4_Table = false;
Question_4_Graph = true;
Question_5_Graph = true;

if Question_1
    disp("Calculating Approximation of B(n, r)");
    trials = 300;
    
    nr_pairs = ...
        {[5, 2],[5, 5],[5, 9],[5,13],[5,16],...
         [10, 2],[10, 7],[10,18],[10,25],[10,35],[10,50],[10,60],[10,75],[10,100],...
         [15, 2],[15, 7],[15,18],[15,25],[15,35],[15,50],[15,60],[15,75],[15,100],...
         [20, 2],[20, 7],[20,18],[20,25],[20,35],[20,50],[20,60],[20,75],[20,100],};
    
    T1 = approxB(nr_pairs, trials);
    writetable(T1, "Q1.txt")
end

if Question_2
    disp("Calculating Approximation of A(n, d)");
    trials = 300;
    nd_pairs = {[5, 1],[5, 2],[5, 3],[5, 4],...
        [10, 2],[10, 4],[10, 6],[10, 8],[10, 9],...
        [15, 3],[15, 5],[15, 7],[15,9],[15,11],[15,13],[15,13]...
        [20, 4],[20, 6],[20, 8],[20, 10],[20,12],[20,14],[20,16],[20,18]};
    
    T2 = approxA(nd_pairs, trials);
    writetable(T2, "Q2.txt")
    
end

if Question_3
    T1 = readtable("Q1.txt");
    T2 = readtable("Q2.txt");
    
    n1 = T1{:, 'n'};
    r1 = T1{:, 'r'};
    d1 = T1{:, 'Max_d'};
    ir1 = log2(r1)./n1;
    ecr1 = (d1-1)./n1;
    
    n2 = T2{:, 'n'};
    r2 = T2{:, 'Max_r'};
    d2 = T2{:, 'd'};
    ir2 = log2(r2)./n2;
    ecr2 = (d2-1)./n2;
    
    figure(1);
    subplot(2,2,1)
    scatter(ecr1, ir1, 'blue')
    axis([0 1 0 1])
    ylabel('Information Rate') 
    xlabel('Error-Control Rate')
    yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    title(["Graph of information rate against Error Control Rate",...
        " using codes generated to approximate of B(n, r)"]);
    grid on
    
    subplot(2,2,2) 
    scatter(ecr2, ir2, 'red')
    axis([0 1 0 1])
    ylabel('Information Rate') 
    xlabel('Error-Control Rate')
    yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    title(["Graph of information rate against Error Control Rate",...
        " using codes generated to approximate of A(n, d)"]);
    grid on
    
    subplot(2,2,[3 4])
    hold on
    scatter(ecr1, ir1, 'blue')
    scatter(ecr2, ir2, 'red')
    axis([0 1 0 1])
    ylabel('Information Rate') 
    xlabel('Error-Control Rate')
    % Here we define B(n,r) to be the maximum of the hamming distance
    % over all codes with length n and size r
    xs = [0 1];    
    
    legend("Approx B(n, r)", "Approx A(n, d)")
    yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    grid on
    hold off
    
    set(gcf,'position',[10,10,1200,900])
    print("Q3 graph.eps", '-depsc');
end

if Question_4_Table
    disp("Calculating Approximation of B(n, r) for linear codes");
    trials = 100;
    nk_pairs = {[5, 1],[5, 2],[5, 3],[5, 4],...
                [10, 2],[10, 4],[10, 6],[10, 8],[10, 10],...
                [15, 2],[15, 4],[15, 6],[15,8],[15,10],[15,12],[15,14],...
                [20, 2],[20, 4],[20,6],[20,8],[20,10],[20,12],[20,14],[20,16]};
            
    
    T3 = approxLinB(nk_pairs, trials);
    writetable(T3, "Q4a.txt");
    
    disp("Calculating Approximation of A(n, d) for linear codes");
    trials = 100;
    nd_pairs = {[5, 1],[5, 2],[5, 3],[5, 4],...
                [10, 2],[10, 4],[10, 6],[10, 8],[10, 9],...
                [15, 3],[15, 5],[15, 7],[15, 9],[15,11],[15,13],...
                [20, 6],[20, 8],[20, 10],[20, 12],[20,14],[20,16],[20,18]};
    
    T4 = approxLinA(nd_pairs, trials);
    writetable(T4, "Q4b.txt");
    
end

if Question_4_Graph
    T3 = readtable("Q4a.txt");
    T4 = readtable("Q4b.txt");
    
    n3 = T3{:, 'n'};
    k3 = T3{:, 'k'};
    d3 = T3{:, 'Max_d'};
    ir3 = k3./n3;
    ecr3 = (d3-1)./n3;
    
    n4 = T4{:, 'n'};
    k4 = T4{:, 'Max_k'};
    d4 = T4{:, 'd'};
    ir4 = k4./n4;
    ecr4 = (d4-1)./n4;
    
    figure(2)
    subplot(2,2,1)
    scatter(ecr3, ir3, 'blue')
    axis([0 1 0 1])
    ylabel('Information Rate') 
    xlabel('Error-Control Rate')
    yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    title(["Graph of Information Rate against Error Control Rate",...
        " using codes generated to approximate $\hat{B}(n, k)$"], 'Interpreter', 'latex');
    grid on
    
    subplot(2,2,2) 
    scatter(ecr4, ir4, 'red')
    axis([0 1 0 1])
    ylabel('Information Rate') 
    xlabel('Error-Control Rate')
    yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    title(["Graph of Information Rate against Error Control Rate",...
        " using codes generated to approximate $\hat{A}(n, d)$"], 'Interpreter', 'latex');
    grid on
    
    subplot(2,2,[3 4])
    hold on
    scatter(ecr3, ir3, 'blue')
    scatter(ecr4, ir4, 'red')
    axis([0 1 0 1])
    ylabel('Information Rate') 
    xlabel('Error-Control Rate')
    % Here we define B(n,r) to be the maximum of the hamming distance
    % over all codes with length n and size r
    xs = [0 1];
    
    
    legend("Approx \hat{B}(n, r)", "Approx \hat{A}(n, d)")
    yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    grid on
    hold off
    
    set(gcf,'position',[10,10,1200,900])
    print("Q4 graph.eps", '-depsc');
end

if Question_5_Graph
    T2 = readtable("Q2.txt");
    T4 = readtable("Q4b.txt");
    
    n2 = T2{:, 'n'};
    r2 = T2{:, 'Max_r'};
    d2 = T2{:, 'd'};
    ir2 = log2(r2)./n2;
    ecr2 = (d2-1)./n2;
    
    n4 = T4{:, 'n'};
    k4 = T4{:, 'Max_k'};
    d4 = T4{:, 'd'};
    ir4 = k4./n4;
    ecr4 = (d4-1)./n4;
    
    ir = [ir2; ir4];
    ecr = [ecr2; ecr4];
    n = [n2; n4];
    
    xs = linspace(0,1, 101);
    ys1 = (1+(xs.*log2(xs)+(1-xs).*log2(1-xs))).*(xs<1/2);
    ys1(1)=1;
    ys1(end)=0;
    ys2 = 1+(0.5*xs.*log2(0.5*xs)+(1-0.5*xs).*log2(1-0.5*xs));
    ys2(1)=1;
    
    figure(3)
    scatter(ecr, ir, [], [(1-n/max(n)) zeros(length(n),2)], 'filled')
    l1 = line(xs, ys1, 'Color', 'blue');
    l2 = line(xs, ys2, 'Color', 'blue', 'LineStyle', '--');
    axis([0 1 0 1])
    ylabel('Information Rate') 
    xlabel('Error-Control Rate')
    yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    title(["Graph of Information Rate against Error Control Rate",...
        " using codes generated to approximate A(n, d) and $\hat{A}(n, d)$"], 'Interpreter', 'latex');
    grid on
    print("Q5 graph.eps", '-depsc');
end

function T = approxB(nr_pairs, trials)
    % pairs will be a list of pairs of (n, r)
    % Trials will be number of attempts
    
    len = length(nr_pairs);
    Ns = zeros(len, 1);
    Rs = zeros(len, 1);
    MinDs = zeros(len, 1);
    for ind = 1:len
        pair = cell2mat(nr_pairs(ind));
        n = pair(1);
        r = pair(2);
        d = 0;
        for k = 1:trials
            d = max(d, minrandcodedist(n,r));
        end
        Ns(ind) = n;
        Rs(ind) = r;
        MinDs(ind) = d;
        disp(join(["Finished calculations for (n, r) = (", num2str(n),... 
            ", ", num2str(r), ")"], ""));
    end

    T = table(Ns, Rs, MinDs);
    T.Properties.VariableNames = {'n', 'r', 'Max_d'};
end

function T = approxA(nd_pairs, trials)
    len = length(nd_pairs);
    Ns = zeros(len, 1);
    Ds = zeros(len, 1);
    Max_r = zeros(len, 1);
    for ind = 1:len
        pair = cell2mat(nd_pairs(ind));
        n = pair(1);
        d = pair(2);
        r = 0;
        for k = 1:trials
            r = max(r, minrandcodespace(n, d, 5));
        end
        Ns(ind) = n;
        Ds(ind) = d;
        Max_r(ind) = r;
        disp(join(["Finished calculations for (n, d) = (", num2str(n),... 
            ", ", num2str(d), ")"], ""));
    end
    
    T = table(Ns, Ds, Max_r);
    T.Properties.VariableNames = {'n', 'd', 'Max_r'};
end

function T = approxLinB(nk_pairs, trials)
    % given n and k find d
    len = length(nk_pairs);
    Ns = zeros(len, 1);
    Ks = zeros(len, 1);
    Max_d = zeros(len, 1);
    for ind = 1:len
        pair = cell2mat(nk_pairs(ind));
        n = pair(1);
        k = pair(2);
        d = 0;
        for p = 1:trials
            d = max(d, minlinrandcodedist(n, k));
        end
        Ns(ind) = n;
        Ks(ind) = k;
        Max_d(ind) = d;
        disp(join(["Finished calculations for (n, k) = (", num2str(n),... 
            ", ", num2str(k), ")"], ""));
    end
    
    T = table(Ns, Ks, Max_d);
    T.Properties.VariableNames = {'n', 'k', 'Max_d'};
end

function T = approxLinA(nd_pairs, trials)
    len = length(nd_pairs);
    Ns = zeros(len, 1);
    Ds = zeros(len, 1);
    Max_k = zeros(len, 1);
    for ind = 1:len
        pair = cell2mat(nd_pairs(ind));
        n = pair(1);
        d = pair(2);
        k = 0;
        for p = 1:trials
            k = max(k, minlinrandcodespace(n, d, 5));
        end
        Ns(ind) = n;
        Ds(ind) = d;
        Max_k(ind) = k;
        disp(join(["Finished calculations for (n, d) = (", num2str(n),... 
            ", ", num2str(d), ")"], ""));
    end
    
    T = table(Ns, Ds, Max_k);
    T.Properties.VariableNames = {'n', 'd', 'Max_k'};
end

function d = dist(x, y)
    % Calculates Hamminng Space between 2 elements
    % of the Hamming Space
    d=0;
    for i = 1:length(x)
        d = d + abs(x(i)-y(i));
    end
end

function md = mindist(code)
    % Make code a cell array {} of vectors
    md = length(cell2mat(code(1)));
    for ix = 1:length(code)
        x = cell2mat(code(ix));
        for iy = ix+1:length(code)
            y = cell2mat(code(iy));
            md = min(md, dist(x,y));
        end
    end
end

function d = minrandcodedist(n, r)
    code = {};
    for a = 1:r
        vec = mod(randi(2,n,1),2);
        code{end+1} = vec;
    end
    d = mindist(code);
end

function bool = checkdist(code, vec, d)
    bool = true;
    for i = 1:length(code)
        vec2 = cell2mat(code(i));
        if dist(vec, vec2) < d
            bool = false;
            break
        end
    end
end

function r = minrandcodespace(n, d, tries)
    code = {zeros(n,1)};
    r = 1;
    successful = true;
    while successful
        i = 0;
        while i <= tries
            if i ~= tries
                % Tries multiple times to find a vector not in the code
                % (i.e not of distance 0) and has distance > d for all
                % elements in code
                vec = mod(randi(2,n,1),2);
                if checkdist(code, vec, d)
                    r = r + 1;
                    code{end+1} = vec;
                    break
                end
            else
                successful = false;
            end
            i = i+1;
        end
    end
end

function md = minlindist(gens)
    % gens will be the generator matrix
    dim = size(gens);
    % codes = zeros(2^dim(1), dim(2));
    md = Inf;
    for i = 1:(2^dim(1)-1)
        s = dec2bin(i);
        s = convertStringsToChars(s);
        vec = zeros(1, dim(1));
        for j = 1:length(s)
            vec(dim(1)+1-j) = str2double(s(length(s)+1-j));
        end
        cw = mod(vec*gens, 2);
        c = dist(cw, zeros(1, dim(2)));
        md = min(md, c);
    end
    % Note that if we get minimum distance of 0
    % then rows of the generator matrix are linearly dependent
    
end

function d = minlinrandcodedist(n, k)
    gen = zeros(k, n);
    for a = 1:k
        vec = mod(randi(2,1,n),2);
        gen(a,:) = vec;
    end
    d = minlindist(gen);
end

function k = minlinrandcodespace(n, d, tries)
    code = zeros(0,n);
    k = 0;
    successful = true;
    while successful
        i = 0;
        while i <= tries
            if i ~= tries
                % Tries multiple times to find a vector not in the code
                % (i.e not of distance 0) and has distance >= d for all
                % elements in code
                vec = mod(randi(2,1,n),2);
                while all(vec == zeros(1,n))
                    % generate a new vec if we get the zero vector
                    vec = mod(randi(2,1,n),2);
                end
                if minlindist([code;vec]) >= d
                    k = k + 1;
                    code = [code;vec];
                    break
                end
            else
                successful = false;
            end
            i = i+1;
        end
    end
end