function outputArg1 = GDe2P(Theta,Phi,y,measurement)
if GDe1Entanglement(0,Theta,Phi,y,measurement)<=0
      outputArg1=0;
      return;
end
Header=1;
func = @(p)GDe1Entanglement(p,Theta,Phi,y,measurement); % function of p alone
options = optimset('FunValCheck','on');
outputArg1 = fzeroTest(func,[0,Header],options);
end
function [b,fval,exitflag,output] = fzeroTest(FunFcnIn,x,options,varargin)
% Initialization
fcount = 0;
iter = 0;
intervaliter = 0;
exitflag = 1;
procedure = ' ';

defaultopt = struct( ...
    'Display', 'notify', ...
    'TolX', eps, ...
    'FunValCheck', 'off', ...
    'OutputFcn', [], ...
    'PlotFcns', []);

% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && strcmpi(FunFcnIn,'defaults')
    b = defaultopt;
    return
end

% initialization
if nargin < 3
   options = []; 
end
% Detect problem structure input
if nargin == 1
    if isa(FunFcnIn,'struct')
        [FunFcnIn,x,options] = separateOptimStruct(FunFcnIn);
    else % Single input and non-structure.
        error('MATLAB:fzero:InputArg',...
            getString(message('MATLAB:optimfun:fzero:InputArg')));
    end
end

if nargin == 0 
    error('MATLAB:fzero:NotEnoughInputs',...
        getString(message('MATLAB:optimfun:fzero:NotEnoughInputs'))); 
end

% Check for non-double inputs
if ~isa(x,'double')
  error('MATLAB:fzero:NonDoubleInput',...
    getString(message('MATLAB:optimfun:fzero:NonDoubleInput')));
end

% Check that options is a struct
if ~isempty(options) && ~isa(options,'struct')
    error('MATLAB:fzero:ArgNotStruct',...
        getString(message('MATLAB:optimfun:commonMessages:ArgNotStruct', 3)));
end

tol = optimget(options,'TolX',defaultopt,'fast');
funValCheck = strcmp(optimget(options,'FunValCheck',defaultopt,'fast'),'on');
printtype = optimget(options,'Display',defaultopt,'fast');
switch printtype
    case {'notify','notify-detailed'}
        trace = 1;
    case {'none', 'off'}
        trace = 0;
    case {'iter','iter-detailed'}
        trace = 3;
    case {'final','final-detailed'}
        trace = 2;
    otherwise
        trace = 1;
end
% Handle the output functions
outputfcn = optimget(options,'OutputFcn',defaultopt,'fast');
if isempty(outputfcn)
    haveoutputfcn = false;
else
    haveoutputfcn = true;
    % Parse OutputFcn which is needed to support cell array syntax for OutputFcn.
    outputfcn = createCellArrayOfFunctions(outputfcn,'OutputFcn');
end

% Handle the plot functions
plotfcns = optimget(options,'PlotFcns',defaultopt,'fast');
if isempty(plotfcns)
    haveplotfcn = false;
else
    haveplotfcn = true;
    % Parse PlotFcns which is needed to support cell array syntax for PlotFcns.
    plotfcns = createCellArrayOfFunctions(plotfcns,'PlotFcns');
end

% Convert to function handle as needed.
[FunFcn,errStruct] = fcnchk(FunFcnIn,length(varargin));
if ~isempty(errStruct)
    error(message(errStruct.identifier));
end
% We know fcnchk succeeded if we got to here
if isa(FunFcn,'inline')      
    if isa(FunFcnIn,'inline')
        Ffcnstr = inputname(1);  % name of inline object such as f where f=inline('x*2');
        if isempty(Ffcnstr)  % inline('sin(x)')  
            Ffcnstr = formula(FunFcn);  % Grab formula, no argument name 
        end
        Ftype = 'inline object';
    else  % not an inline originally (character array expression).
        Ffcnstr = char(FunFcnIn);  % get the character array expression
        Ftype = 'expression';
    end
elseif isa(FunFcn,'function_handle') % function handle
    Ffcnstr = func2str(FunFcn);  % get the name passed in
    Ftype = 'function_handle';
else  % Not converted, must be m-file or builtin
    Ffcnstr = char(FunFcnIn);  % get the name passed in
    Ftype = 'function';
end

% Add a wrapper function to check for Inf/NaN/complex values
if funValCheck
    % Add a wrapper function, CHECKFUN, to check for NaN/complex values without
    % having to change the calls that look like this:
    % f = funfcn(x,varargin{:});
    % x is the first argument to CHECKFUN, then the user's function,
    % then the elements of varargin. To accomplish this we need to add the 
    % user's function to the beginning of varargin, and change funfcn to be
    % CHECKFUN.
    varargin = [{FunFcn}, varargin];
    FunFcn = @checkfun;
end

% Initialize the output and plot functions.
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,[],'init',fcount,iter,intervaliter, ...
        [],procedure,[],[],[],[],varargin{:});
    if stop
        [b,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
        if  trace > 0
            disp(output.message)
        end
        return;
    end
end

if  ~all(isfinite(x))
    error('MATLAB:fzero:Arg2NotFinite',...
        getString(message('MATLAB:optimfun:fzero:Arg2NotFinite')));
end

% Interval input
if (numel(x) == 2) 
    if trace > 2
        disp(' ') %Initial blank line
    end
    a = x(1); savea=a;
    b = x(2); saveb=b;
    % Put first feval in try catch
    try
        fa = FunFcn(a,varargin{:});
    catch ME
        if ~isempty(Ffcnstr)
            error('MATLAB:fzero:InvalidFunctionSupplied',...
                getString(message('MATLAB:optimfun:fzero:InvalidFunctionSupplied',sprintf('%s ==> %s',Ftype,Ffcnstr),ME.message)));
        else
            error('MATLAB:fzero:InvalidFunctionSupplied',...
                getString(message('MATLAB:optimfun:fzero:InvalidFunctionSupplied',Ftype,ME.message)));
        end
        
    end
    
    fb = FunFcn(b,varargin{:});
    if any(~isfinite([fa fb])) || any(~isreal([fa fb]))
        error('MATLAB:fzero:ValuesAtEndPtsComplexOrNotFinite',...
            getString(message('MATLAB:optimfun:fzero:ValuesAtEndPtsComplexOrNotFinite')));
    end
    fcount = fcount + 2;
    savefa = fa; savefb = fb;
    
    if ( fa == 0 )
        b = a;
        msg = getString(message('MATLAB:optimfun:fzero:ZeroFindTerminated'));
        if trace > 1
            disp(msg)
        end
        output.intervaliterations = intervaliter;
        output.iterations = iter;
        output.funcCount = fcount;
        output.algorithm = 'bisection, interpolation';
        output.message = msg;
        fval = fa;
        return
    elseif ( fb == 0)
        % b = b;
        msg = getString(message('MATLAB:optimfun:fzero:ZeroFindTerminated'));
        if trace > 1
            disp(msg)
        end
        output.intervaliterations = intervaliter;
        output.iterations = iter;
        output.funcCount = fcount;
        output.algorithm = 'bisection, interpolation';
        output.message = msg;
        fval = fb;
        return
    elseif (fa > 0) == (fb > 0)
        error('MATLAB:fzero:ValuesAtEndPtsSameSign',...
            getString(message('MATLAB:optimfun:fzero:ValuesAtEndPtsSameSign')));
    end
    
    % Starting guess scalar input
elseif (numel(x) == 1)
    if trace > 2 
        disp(' ')
        fprintf(getString(message('MATLAB:optimfun:fzero:SearchForIntervalContainingSignChange',sprintf('%g',x))));
        header = ' Func-count    a          f(a)             b          f(b)        Procedure';
    end
    % Put first feval in try catch
    try
        fx = FunFcn(x,varargin{:});
    catch ME
        if ~isempty(Ffcnstr)
            error('MATLAB:fzero:InvalidFunctionSupplied',...
                getString(message('MATLAB:optimfun:fzero:InvalidFunctionSupplied', sprintf('%s ==> %s',Ftype,Ffcnstr),ME.message)));
        else
            error('MATLAB:fzero:InvalidFunctionSupplied',...
                getString(message('MATLAB:optimfun:fzero:InvalidFunctionSupplied',Ftype,ME.message)));
        end   
    end
    fcount = fcount + 1;  
    if fx == 0
        b = x;
        msg = getString(message('MATLAB:optimfun:fzero:ZeroFindTerminated'));
        if trace > 1
            disp(msg)
        end
        output.intervaliterations = intervaliter;
        output.iterations = iter;
        output.funcCount = fcount;
        output.algorithm = 'bisection, interpolation';
        output.message = msg;
        fval = fx;
        return
    elseif ~isfinite(fx) || ~isreal(fx)
        error('MATLAB:fzero:ValueAtInitGuessComplexOrNotFinite',...
            getString(message('MATLAB:optimfun:fzero:ValueAtInitGuessComplexOrNotFinite')));
    end
    
    if x ~= 0
        dx = x/50;
    else 
        dx = 1/50;
    end
    
    % Find change of sign.
    twosqrt = sqrt(2); 
    a = x; fa = fx; b = x; fb = fx;
    
    if trace > 2
        disp(header)
        procedure='initial interval';
        fprintf('%5.0f   %13.6g %13.6g %13.6g %13.6g   %s\n',fcount,a,fa,b,fb, procedure);
    end
    % OutputFcn and PlotFcns call
    if haveoutputfcn || haveplotfcn
        [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,x,'iter',fcount,iter,intervaliter, ...
            fx,procedure,a,fa,b,fb,varargin{:}); % a and b are x to start
        if stop
            [b,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
            if  trace > 0
                disp(output.message)
            end
            return;
        end
    end

    while (fa > 0) == (fb > 0)
        intervaliter = intervaliter + 1;
        dx = twosqrt*dx;
        a = x - dx;  fa = FunFcn(a,varargin{:});
        fcount = fcount + 1;
        if ~isfinite(fa) || ~isreal(fa) || ~isfinite(a)
            [exitflag,msg] = disperr(a,fa,trace);
            b = NaN; fval = NaN;
            output.intervaliter = intervaliter;
            output.iterations = iter;
            output.funcCount = fcount;
            output.algorithm = 'bisection, interpolation';
            output.message = msg;
            return
        end

        if (fa > 0) ~= (fb > 0) % check for different sign
            % Before we exit the while loop, print out the latest interval
            if trace > 2
                procedure='search';
                fprintf('%5.0f   %13.6g %13.6g %13.6g %13.6g   %s\n',fcount,a,fa,b,fb, procedure);
            end
            % OutputFcn and PlotFcns call
            if haveoutputfcn || haveplotfcn
                [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,x,'iter',fcount,iter,intervaliter, ...
                    fx,procedure,a,fa,b,fb,varargin{:});
                if stop
                    [b,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
                    if  trace > 0
                        disp(output.message)
                    end
                    return;
                end
            end
            break
        end
        
        b = x + dx;  fb = FunFcn(b,varargin{:});
        if ~isfinite(fb) || ~isreal(fb) || ~isfinite(b)
            [exitflag,msg] = disperr(b,fb,trace);
            b = NaN; fval = NaN;
            output.intervaliter = intervaliter;
            output.iterations = iter;
            output.funcCount = fcount;
            output.algorithm = 'bisection, interpolation';
            output.message = msg;
            return
        end
        fcount = fcount + 1;        
        if trace > 2
            procedure='search';
            fprintf('%5.0f   %13.6g %13.6g %13.6g %13.6g   %s\n',fcount,a,fa,b,fb, procedure);
        end
        % OutputFcn and PlotFcns call
        if haveoutputfcn || haveplotfcn
            [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,x,'iter',fcount,iter,intervaliter, ...
                fx,procedure,a,fa,b,fb,varargin{:});
            if stop
                [b,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
                if  trace > 0
                    disp(output.message)
                end
                return;
            end
        end
    end % while
    
    if trace > 2
        disp(' ')
        fprintf(getString(message('MATLAB:optimfun:fzero:SearchForZeroInInterval',sprintf('%g',a),sprintf('%g',b))));
    end
    savea = a; savefa = fa; saveb = b; savefb = fb;
else
    error('MATLAB:fzero:LengthArg2',...
        getString(message('MATLAB:optimfun:fzero:LengthArg2')));
end % if (numel(x) == 2)

fc = fb;
procedure = 'initial';
header2 = ' Func-count    x          f(x)             Procedure';
if trace > 2
    disp(header2)
end
% Main loop, exit from middle of the loop
while fb ~= 0 && a ~= b
    % Insure that b is the best result so far, a is the previous
    % value of b, and c is on the opposite side of the zero from b.
    if (fb > 0) == (fc > 0)
        c = a;  fc = fa;
        d = b - a;  e = d;
    end
    if abs(fc) < abs(fb)
        a = b;    b = c;    c = a;
        fa = fb;  fb = fc;  fc = fa;
    end
    
    % Convergence test and possible exit
    m = 0.5*(c - b);
    toler = 2.0*tol*max(abs(b),1.0);
    if (abs(m) <= toler) || (fb == 0.0) 
        break
    end
    if trace > 2
        fprintf('%5.0f   %13.6g %13.6g        %s\n',fcount, b, fb, procedure);
    end
    % OutputFcn and PlotFcns call
    if haveoutputfcn || haveplotfcn
        [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,b,'iter',fcount,iter,intervaliter, ...
            fb,procedure,savea,savefa,saveb,savefb,varargin{:});
        if stop
            [b,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
            if  trace > 0
                disp(output.message)
            end
            return;
        end
    end
    
    % Choose bisection or interpolation
    if (abs(e) < toler) || (abs(fa) <= abs(fb))
        % Bisection
        d = m;  e = m;
        procedure='bisection';
    else
        % Interpolation
        s = fb/fa;
        if (a == c)
            % Linear interpolation
            p = 2.0*m*s;
            q = 1.0 - s;
        else
            % Inverse quadratic interpolation
            q = fa/fc;
            r = fb/fc;
            p = s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
            q = (q - 1.0)*(r - 1.0)*(s - 1.0);
        end
        if p > 0
            q = -q;
        else
            p = -p;
        end
        % Is interpolated point acceptable
        if (2.0*p < 3.0*m*q - abs(toler*q)) && (p < abs(0.5*e*q))
            e = d;  d = p/q;
            procedure='interpolation';
        else
            d = m;  e = m;
            procedure='bisection';
        end
    end % Interpolation
    
    % Next point
    a = b;
    fa = fb;
    if abs(d) > toler
        b = b + d;
    elseif b > c
        b = b - toler;
    else
        b = b + toler;
    end
    fb = FunFcn(b,varargin{:});
    fcount = fcount + 1;
    iter = iter + 1;
end % Main loop

fval = fb; % b is the best value

% Output last chosen b
if trace > 2
    fprintf('%5.0f   %13.6g %13.6g        %s\n',fcount, b, fb, procedure);
end

% OutputFcn and PlotFcns call
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,b,'iter',fcount,iter,intervaliter, ...
        fb,procedure,savea,savefa,saveb,savefb,varargin{:});
    if stop
        [b,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
        if  trace > 0
            disp(output.message)
        end
        return;
    end
end

output.intervaliterations = intervaliter;
output.iterations = iter;
output.funcCount = fcount;
output.algorithm = 'bisection, interpolation';

if abs(fval) <= max(abs(savefa),abs(savefb))
    msg = sprintf(getString(message('MATLAB:optimfun:fzero:ZeroFoundInInterval',sprintf('%g',savea),sprintf('%g',saveb))));
else
    exitflag = -5; 
    msg = sprintf(...
        getString(message('MATLAB:optimfun:fzero:CurrentPointXNearSingularPoint',...
        sprintf('%g',savea),sprintf('%g',saveb))));
end
if trace > 1
    disp(' ')
    disp(msg)
end
output.message = msg;
% Outputfcn and PlotFcns call
if haveoutputfcn || haveplotfcn
    callOutputAndPlotFcns(outputfcn,plotfcns,b,'done',fcount,iter,intervaliter,fval,procedure,savea,savefa,saveb,savefb,varargin{:});
end

end
%------------------------------------------------------------------
function [exitflag,msg] = disperr(y, fy, trace)
%DISPERR Display an appropriate error message when FY is Inf, 
%   NaN, or complex.  Assumes Y is the value and FY is the function 
%   value at Y. If FY is neither Inf, NaN, or complex, it generates 
%   an error message.

if ~isfinite(fy)  % NaN or Inf detected
    exitflag = -3;
    msg = ...
        getString(message('MATLAB:optimfun:fzero:AbortingSearchForIntervalNaNInfEncountered',sprintf('%g',y),sprintf('%g',fy)));
    if trace > 0
        disp(msg)
    end
elseif ~isreal(fy) % Complex value detected
    exitflag = -4;
    msg = ...
        getString(message('MATLAB:optimfun:fzero:AbortingSearchForIntervalComplexEncountered',sprintf('%g',y),num2str(fy)));
    if trace > 0
        disp(msg)        
    end
elseif ~isfinite(y) % Inf detected in bracketting stage
    exitflag = -6;
    msg = ...
        getString(message('MATLAB:optimfun:fzero:ExitingNoSignChange'));
    if trace > 0
        disp(msg)        
    end
end
end
%--------------------------------------------------------------------------
function [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,x,state,fcount,iter,intervaliter,  ...
    f,procedure,a,fvala,b,fvalb,varargin)
% CALLOUTPUTANDPLOTFCNS assigns values to the struct OptimValues and then calls the
% outputfcn/plotfcns.  
%
% state - can have the values 'init','iter', or 'done'. 
% We do not handle the case 'interrupt' because we do not want to update
% xOutputfcn or optimValues (since the values could be inconsistent) before calling
% the outputfcn; in that case the outputfcn/plotfcns are called directly rather than
% calling it inside callOutputAndPlotFcns.

% For the 'done' state we do not check the value of 'stop' because the
% optimization is already done.
optimValues.funccount = fcount;
optimValues.iteration = iter;
optimValues.intervaliteration = intervaliter;
optimValues.fval = f;
optimValues.procedure = procedure;
optimValues.intervala = a;
optimValues.fvala = fvala;
optimValues.intervalb = b;
optimValues.fvalb = fvalb;

xOutputfcn = x;  % set xOutputfcn to be x
stop = false;
% Call output functions
if ~isempty(outputfcn)
    switch state
        case {'iter','init'}
            stop = callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:});
    end
end
% Call plot functions
if ~isempty(plotfcns)
    switch state
        case {'iter','init'}
            stop = callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:});
    end
end
end
%--------------------------------------------------------------------------
function [b,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues)
% CLEANUPINTERRUPT updates or sets all the output arguments of FMINBND when the optimization 
% is interrupted.

% Call plot function driver to finalize the plot function figure window. If
% no plot functions have been specified or the plot function figure no
% longer exists, this call just returns.
callAllOptimPlotFcns('cleanuponstopsignal');

b = xOutputfcn;
fval = optimValues.fval;
exitflag = -1; 
output.intervaliterations = optimValues.intervaliteration;
output.iterations = optimValues.iteration;
output.funcCount = optimValues.funccount;
output.algorithm = 'bisection, interpolation';
output.message = getString(message('MATLAB:optimfun:fzero:OptimizationTerminatedPrematurelyByUser'));
end
%--------------------------------------------------------------------------
function f = checkfun(x,userfcn,varargin)
% CHECKFUN checks for complex or NaN results from userfcn.

f = userfcn(x,varargin{:});
% Note: we do not check for Inf as FZERO handles it naturally.  ???
if isnan(f)
    error('MATLAB:fzero:checkfun:NaNFval',...
        getString(message('MATLAB:optimfun:fzero:checkfun:NaNFval', localChar( userfcn ), sprintf( '%g', x ))));  
elseif ~isreal(f)
    error('MATLAB:fzero:checkfun:ComplexFval',...
        getString(message('MATLAB:optimfun:fzero:checkfun:ComplexFval', localChar( userfcn ), sprintf( '%g', x ))));  
end
end
%--------------------------------------------------------------------------
function strfcn = localChar(fcn)
% Convert the fcn to a character array for printing

if ischar(fcn)
    strfcn = fcn;
elseif isa(fcn,'inline')
    strfcn = char(fcn);
elseif isa(fcn,'function_handle')
    strfcn = func2str(fcn);
else
    try
        strfcn = char(fcn);
    catch ME
        strfcn = getString(message('MATLAB:optimfun:fzero:NameNotPrintable'));
    end
end
end