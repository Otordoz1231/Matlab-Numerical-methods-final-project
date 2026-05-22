function NumericalMethods()


fig = figure('Name','Numerical Methods PRO', ...
             'NumberTitle','off', ...
             'Position',[100 80 1150 720], ...
             'Resize','on');

tg = uitabgroup(fig,'Position',[0 0 1 1]);
t1 = uitab(tg,'Title','Root Finding');
t2 = uitab(tg,'Title','Matrix');

buildRootFindingTab(t1);
buildMatrixTab(t2);
end


% =========================================================
%  ROOT FINDING TAB
% =========================================================
function buildRootFindingTab(parent)

% --- equation entry
uicontrol(parent,'Style','text','String','Equation f(x)', ...
    'Position',[20 665 120 20],'HorizontalAlignment','left');

hEq = uicontrol(parent,'Style','edit','String','x.^3-x-2', ...
    'Position',[140 663 400 24]);

% --- method selector
uicontrol(parent,'Style','text','String','Method:', ...
    'Position',[560 665 60 20],'HorizontalAlignment','left');

methods = {'Incremental','Bisection','False Position','Newton','Secant'};
hMethod = uicontrol(parent,'Style','popupmenu','String',methods, ...
    'Value',2, ...
    'Position',[620 663 160 24], ...
    'Callback',@(~,~) updateInputs());

% --- Tolerance input
uicontrol(parent,'Style','text','String','Tolerance:', ...
    'Position',[800 665 70 20],'HorizontalAlignment','left');
hTol = uicontrol(parent,'Style','edit','String','0.001', ...
    'Position',[875 663 80 24]);

% --- Max Iterations input
uicontrol(parent,'Style','text','String','Max Iter:', ...
    'Position',[970 665 65 20],'HorizontalAlignment','left');
hMaxIter = uicontrol(parent,'Style','edit','String','50', ...
    'Position',[1040 663 80 24]);

% --- column header label (blue, bold - matches Python)
hHeader = uicontrol(parent,'Style','text','String','', ...
    'Position',[20 638 1100 20], ...
    'ForegroundColor','blue','FontWeight','bold', ...
    'HorizontalAlignment','left');

% --- dynamic input panel
hInputPanel = uipanel(parent,'Title','Inputs', ...
    'Position',[0.01 0.80 0.55 0.10]);

inputHandles = struct();

% --- results table
hTable = uitable(parent, ...
    'Position',[20 65 1110 490], ...
    'RowName',{}, ...
    'ColumnName',{});

% --- Solve button (green)
uicontrol(parent,'Style','pushbutton','String','Solve', ...
    'Position',[20 20 100 35], ...
    'BackgroundColor',[0 0.6 0], ...
    'ForegroundColor','w','FontWeight','bold', ...
    'Callback',@(~,~) solveCB());

% --- Plot Graph button (blue)
uicontrol(parent,'Style','pushbutton','String','Plot Graph', ...
    'Position',[140 20 120 35], ...
    'BackgroundColor',[0 0 0.8], ...
    'ForegroundColor','w','FontWeight','bold', ...
    'Callback',@(~,~) graphCB());

updateInputs();

    % ==========================================================
    %  updateInputs  - matches Python's update() exactly
    %  Field names and headers are 1:1 with the Python code
    % ==========================================================
    function updateInputs()
        delete(get(hInputPanel,'Children'));
        inputHandles = struct();

        m = methods{get(hMethod,'Value')};

        switch m
            case 'Bisection'
                hHeader.String = ...
                  'Iteration | XL | XR | XU | f(XL) | f(XR) | Ea | Product | Remarks';
                fields    = {'XL','XU'};
                labels    = {'XL','XU'};

            case 'False Position'
                hHeader.String = ...
                  'Iteration | XL | XU | XR | Ea | f(XL) | f(XU) | f(XR) | Product';
                fields    = {'XL','XU'};
                labels    = {'XL','XU'};

            case 'Secant'
                hHeader.String = ...
                  'Iteration | Xi-1 | Xi | Xi+1 | Ea';
                fields    = {'Xi_prev','Xi'};
                labels    = {'Xi-1','Xi'};

            case 'Newton'
                hHeader.String = 'Iteration | Xi | Ea | f(x) | f''(x)';
                fields    = {'Xi'};
                labels    = {'Xi'};

            case 'Incremental'
                hHeader.String = ...
                  'Iteration | XL | Step | XU | f(XL) | f(XU) | Product | Remarks';
                fields    = {'XL','Step'};
                labels    = {'XL','Step'};
        end

        for k = 1:numel(fields)
            xpos = 20 + (k-1)*200;
            uicontrol(hInputPanel,'Style','text','String',labels{k}, ...
                'Position',[xpos 32 60 18],'HorizontalAlignment','center');
            he = uicontrol(hInputPanel,'Style','edit','String','', ...
                'Position',[xpos 10 80 22]);
            inputHandles.(fields{k}) = he;
        end
    end

    % ==========================================================
    %  solveCB  - table columns match Python's cols tuples exactly
    % ==========================================================
    function solveCB()
        try
            eqStr   = get(hEq,'String');
            m       = methods{get(hMethod,'Value')};
            tol     = getTol();
            maxIter = getMaxIter();

            switch m

                case 'Bisection'
                    xl = getVal('XL'); xu = getVal('XU');
                    [rootv, T] = do_bisection(eqStr, xl, xu, tol, maxIter);
                    cols = {'i','XL','XR','XU','f(XL)','f(XR)','Ea','Product','Remarks'};

                case 'False Position'
                    xl = getVal('XL'); xu = getVal('XU');
                    [rootv, T] = do_false_position(eqStr, xl, xu, tol, maxIter);
                    cols = {'i','XL','XU','XR','Ea','f(XL)','f(XU)','f(XR)','Product','Remarks'};

                case 'Newton'
                    x0 = getVal('Xi');
                    [rootv, T] = do_newton(eqStr, x0, tol, maxIter);
                    cols = {'i','Xi','Ea','f(x)','f''(x)'};

                case 'Secant'
                    x0 = getVal('Xi_prev');
                    x1 = getVal('Xi');
                    [rootv, T] = do_secant(eqStr, x0, x1, tol, maxIter);
                    cols = {'i','Xi-1','Xi','Xi+1','Ea','f(Xi-1)','f(Xi)','f(Xi+1)'};

                case 'Incremental'
                    xl   = getVal('XL');
                    step = getVal('Step');
                    [rootv, T] = do_incremental(eqStr, xl, step, tol, maxIter);
                    cols = {'i','XL','Step','XU','f(XL)','f(XU)','Product','Remarks'};
            end

            set(hTable,'ColumnName',cols,'Data',T);
            msgbox(sprintf('Root  %.3f', rootv),'Result');

        catch ME
            errordlg(ME.message,'Error');
        end
    end

    % ==========================================================
    %  graphCB  - mirrors Python graph() function exactly
    % ==========================================================
    function graphCB()
        try
            eqStr   = get(hEq,'String');
            m       = methods{get(hMethod,'Value')};
            tol     = getTol();
            maxIter = getMaxIter();

            switch m

                case 'Bisection'
                    xl = getVal('XL'); xu = getVal('XU');
                    if xl > xu, [xl,xu] = deal(xu,xl); end
                    [rootv,~] = do_bisection(eqStr, xl, xu, tol, maxIter);
                    pts  = {'XL',xl,'blue'; 'XU',xu,'green'; 'XR',rootv,'red'};
                    xmin = xl-2; xmax = xu+2;

                case 'False Position'
                    xl = getVal('XL'); xu = getVal('XU');
                    if xl > xu, [xl,xu] = deal(xu,xl); end
                    [rootv,~] = do_false_position(eqStr, xl, xu, tol, maxIter);
                    pts  = {'XL',xl,'blue'; 'XU',xu,'green'; 'XR',rootv,'red'};
                    xmin = xl-2; xmax = xu+2;

                case 'Newton'
                    x0 = getVal('Xi');
                    [rootv,~] = do_newton(eqStr, x0, tol, maxIter);
                    pts  = {'Xi',x0,'blue'; 'Root',rootv,'red'};
                    xmin = x0-5; xmax = x0+5;

                case 'Secant'
                    x0 = getVal('Xi_prev');
                    x1 = getVal('Xi');
                    [rootv,~] = do_secant(eqStr, x0, x1, tol, maxIter);
                    pts  = {'Xi-1',x0,'blue'; 'Xi',x1,'green'; 'Root',rootv,'red'};
                    xmin = min(x0,x1)-2; xmax = max(x0,x1)+2;

                case 'Incremental'
                    xl   = getVal('XL');
                    step = getVal('Step');
                    [rr,~] = do_incremental(eqStr, xl, step, tol, maxIter);
                    if isempty(rr)
                        errordlg('No bracket found within max iterations.','Graph Error');
                        return;
                    end
                    rootv = mean(rr);
                    pts   = {'XL',rr(1),'blue'; 'XU',rr(2),'green'; 'Root',rootv,'red'};
                    xmin  = xl-2; xmax = xl+step*10;
            end

            xv = linspace(xmin, xmax, 1000);
            yv = arrayfun(@(xval) feval_safe(eqStr, xval), xv);

            figure('Name',[m ' Method Graph']);
            plot(xv, yv, 'LineWidth', 2, 'DisplayName','f(x)');
            hold on;
            yline(0, 'k');

            for k = 1:size(pts,1)
                pname  = pts{k,1};
                px     = pts{k,2};
                pcolor = pts{k,3};
                py     = feval_safe(eqStr, px);
                scatter(px, py, 120, pcolor, 'filled');
                text(px, py, sprintf(' %s\n(%.3f, %.3f)', pname, px, py), ...
                    'FontSize',9);
            end

            grid on;
            title([m ' Method Graph']);
            xlabel('x'); ylabel('f(x)');
            legend('f(x)');

        catch ME
            errordlg(ME.message,'Graph Error');
        end
    end

    % --- helper: read a numeric entry by field key
    function v = getVal(field)
        h = inputHandles.(field);
        v = str2double(get(h,'String'));
        if isnan(v)
            error('Invalid input for "%s". Please enter a number.', field);
        end
    end

    % --- helper: read tolerance
    function v = getTol()
        v = str2double(get(hTol,'String'));
        if isnan(v) || v <= 0
            error('Tolerance must be a positive number (e.g. 0.001).');
        end
    end

    % --- helper: read max iterations
    function v = getMaxIter()
        v = round(str2double(get(hMaxIter,'String')));
        if isnan(v) || v < 1
            error('Max Iterations must be a positive integer (e.g. 50).');
        end
    end

end  % buildRootFindingTab


% =========================================================
%  MATRIX TAB
% =========================================================
function buildMatrixTab(parent)

uicontrol(parent,'Style','text', ...
    'String','Matrix A (example: 1 2;3 4)', ...
    'Position',[20 668 400 20],'HorizontalAlignment','left');
hA = uicontrol(parent,'Style','edit','String','1 2;3 4', ...
    'Position',[20 644 450 26]);

uicontrol(parent,'Style','text','String','Matrix B', ...
    'Position',[20 618 100 20],'HorizontalAlignment','left');
hB = uicontrol(parent,'Style','edit','String','5 6;7 8', ...
    'Position',[20 594 450 26]);

hOut = uicontrol(parent,'Style','listbox','String',{}, ...
    'Position',[20 60 1100 520], ...
    'FontName','Courier','FontSize',10, ...
    'BackgroundColor','w');

% Operations matching Python ops dict exactly
ops = { 'Add',         @(a,b) a+b; ...
        'Multiply',    @(a,b) a*b; ...
        'Transpose',   @(a,b) a'; ...
        'Determinant', @(a,b) det(a); ...
        'Inverse',     @(a,b) inv(a); ...
        'Adjoint',     @(a,b) inv(a)*det(a); ...
        'Power',       @(a,b) a^2; ...
        'Solve Ax=b',  @(a,b) a\b };

bx = 20;
for k = 1:size(ops,1)
    uicontrol(parent,'Style','pushbutton','String',ops{k,1}, ...
        'Position',[bx 20 120 32], ...
        'Callback',@(~,~) doOp(ops{k,2}));
    bx = bx + 130;
end

    function doOp(func)
        try
            A = parseMatrix(get(hA,'String'));
            B = parseMatrix(get(hB,'String'));
            result = func(A, B);
            lines  = formatMatrix(result);
            set(hOut,'String',lines);
        catch ME
            errordlg(ME.message,'Matrix Error');
        end
    end

    function M = parseMatrix(txt)
        rows = strsplit(strtrim(txt), ';');
        M = cellfun(@(r) str2num(strtrim(r)), rows, 'UniformOutput', false); %#ok<ST2NM>
        M = cell2mat(M');
        if isempty(M)
            error('Format:  1 2;3 4');
        end
    end

    function lines = formatMatrix(x)
        if isscalar(x)
            lines = {num2str(x, '%.6g')};
        elseif isvector(x)
            lines = cellstr(num2str(x(:), '%.6g'));
        else
            lines = cell(size(x,1), 1);
            for r = 1:size(x,1)
                lines{r} = num2str(x(r,:), '%12.6g');
            end
        end
    end

end  % buildMatrixTab


% =========================================================
%  HELPER: safe equation evaluator  (supports sin, cos, exp, log, etc.)
% =========================================================
function y = feval_safe(eq, x)
    y = eval(eq);
end

% =========================================================
%  HELPER: numerical derivative (central difference, h=1e-5)
%  Mirrors Python:  (f(x+h) - f(x-h)) / (2h)
% =========================================================
function dy = numDeriv(eq, x)
    h  = 1e-5;
    dy = (feval_safe(eq, x+h) - feval_safe(eq, x-h)) / (2*h);
end

% =========================================================
%  HELPER: stopping criterion
%  Mirrors Python stop():  abs(val)<tol  OR  (ea!=0 AND ea<tol)
% =========================================================
function s = shouldStop(ea, val, tol)
    s = abs(val) < tol || (ea ~= 0 && ea < tol);
end

% =========================================================
%  HELPER: round to 3 decimal places  (mirrors Python r())
% =========================================================
function v = r3(v)
    v = round(v, 3);
end


% =========================================================
%  BISECTION
%  Table columns: i, XL, XR, XU, f(XL), f(XR), Ea, Product, Remarks
% =========================================================
function [xr, T] = do_bisection(eq, xl, xu, tol, maxIter)
    if xl > xu, [xl, xu] = deal(xu, xl); end
    if feval_safe(eq,xl) * feval_safe(eq,xu) > 0
        error('f(XL) and f(XU) must have opposite signs');
    end

    T      = {};
    xr_old = [];
    xr     = xl;

    for i = 1:maxIter
        xr  = (xl + xu) / 2;
        fxl = feval_safe(eq, xl);
        fxr = feval_safe(eq, xr);

        if isempty(xr_old)
            ea = 0;
        else
            ea = abs((xr - xr_old) / xr) * 100;
        end

        prod_val = fxl * fxr;

        if prod_val < 0
            rem = '<0 Revert back to XL';
        else
            rem = '>0 Go to next interval';
        end

        T(end+1, :) = {i, r3(xl), r3(xr), r3(xu), ...
                       r3(fxl), r3(fxr), r3(ea), r3(prod_val), rem}; %#ok<AGROW>

        if shouldStop(ea, fxr, tol), break; end

        if prod_val < 0
            xu = xr;
        else
            xl = xr;
        end

        xr_old = xr;
    end
end


% =========================================================
%  FALSE POSITION
%  Table columns: i, XL, XU, XR, Ea, f(XL), f(XU), f(XR), Product, Remarks
% =========================================================
function [xr, T] = do_false_position(eq, xl, xu, tol, maxIter)
    if xl > xu, [xl, xu] = deal(xu, xl); end
    if feval_safe(eq,xl) * feval_safe(eq,xu) > 0
        error('f(XL) and f(XU) must have opposite signs');
    end

    T      = {};
    xr_old = [];
    xr     = xl;

    for i = 1:maxIter
        fxl = feval_safe(eq, xl);
        fxu = feval_safe(eq, xu);
        xr  = xu - (fxu * (xl - xu)) / (fxl - fxu);
        fxr = feval_safe(eq, xr);

        if isempty(xr_old)
            ea = 0;
        else
            ea = abs((xr - xr_old) / xr) * 100;
        end

        prod_val = fxl * fxr;

        if prod_val < 0
            rem = '<0 Revert back to XL';
        else
            rem = '>0 Go to next interval';
        end

        T(end+1, :) = {i, r3(xl), r3(xu), r3(xr), r3(ea), ...
                       r3(fxl), r3(fxu), r3(fxr), r3(prod_val), rem}; %#ok<AGROW>

        if shouldStop(ea, fxr, tol), break; end

        if prod_val < 0
            xu = xr;
        else
            xl = xr;
        end

        xr_old = xr;
    end
end


% =========================================================
%  NEWTON-RAPHSON
%  Table columns: i, Xi, Ea, f(x), f'(x)
%  NOTE: Python appends row BEFORE updating x0,
%        and stops on shouldStop(ea, fx)  [checks fx, not f(x1)]
% =========================================================
function [x1, T] = do_newton(eq, x0, tol, maxIter)
    T  = {};
    x1 = x0;

    for i = 1:maxIter
        fx  = feval_safe(eq, x0);
        dfx = numDeriv(eq, x0);

        if dfx == 0, break; end

        x1 = x0 - fx / dfx;

        ea = abs((x1 - x0) / x1) * 100;

        % Row: i, x0 (current Xi), Ea, f(x0), f'(x0)
        T(end+1, :) = {i, r3(x0), r3(ea), r3(fx), r3(dfx)}; %#ok<AGROW>

        if shouldStop(ea, fx, tol), break; end

        x0 = x1;
    end
end


% =========================================================
%  SECANT
%  Table columns: i, Xi-1, Xi, Xi+1, Ea, f(Xi-1), f(Xi), f(Xi+1)
%  NOTE: Python computes ea = abs((x2-x1)/x2)*100  (no zero-guard on first iter
%        because x2 is always computed before ea)
% =========================================================
function [x2, T] = do_secant(eq, x0, x1, tol, maxIter)
    T  = {};
    x2 = x1;   % safe default

    for i = 1:maxIter
        f0 = feval_safe(eq, x0);
        f1 = feval_safe(eq, x1);

        if f0 == f1, break; end

        x2 = x1 - f1 * (x0 - x1) / (f0 - f1);
        f2 = feval_safe(eq, x2);

        ea = abs((x2 - x1) / x2) * 100;

        % Row: i, x0(Xi-1), x1(Xi), x2(Xi+1), Ea, f0, f1, f2
        T(end+1, :) = {i, r3(x0), r3(x1), r3(x2), r3(ea), ...
                       r3(f0), r3(f1), r3(f2)}; %#ok<AGROW>

        if shouldStop(ea, f2, tol), break; end

        x0 = x1;
        x1 = x2;
    end
end


% =========================================================
%  INCREMENTAL SEARCH
%  Table columns: i, XL, Step, XU, f(XL), f(XU), Product, Remarks
%  Returns bracket [x0, x1] when sign change found, else []
%  rootv returned = x1 (last XU before/at bracket) for msgbox
% =========================================================
function [rootv, T] = do_incremental(eq, x0, step, tol, maxIter)
    T     = {};
    rootv = [];      % will be set to x1 of bracket row

    for i = 1:maxIter
        x1 = x0 + step;
        f0 = feval_safe(eq, x0);
        f1 = feval_safe(eq, x1);

        prod_val = f0 * f1;

        if prod_val < 0
            rem = '<0 Revert back to XL';
        else
            rem = '>0 Go to next interval';
        end

        T(end+1, :) = {i, r3(x0), r3(step), r3(x1), ...
                       r3(f0), r3(f1), r3(prod_val), rem}; %#ok<AGROW>

        if abs(f1) < tol || prod_val < 0
            rootv = [x0, x1];   % return bracket pair
            return;
        end

        x0 = x1;
    end
    % rootv stays [] if no bracket found
end
