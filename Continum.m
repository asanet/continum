classdef Continum < handle
    %% Continum: Software for Numerical Continuation and Dynamic Analysis
    % Requires: --- Symbolic math toolbox
    %           --- Optimization toolbox
    %           --- Functions in (provided) 'Aux' folder must be in path
    
    properties (Access = public)
        model, ssmodel, massMatrix, vIndexVar
        valueX, valueP, lowerX, lowerP, upperX, upperP, fixedP, tagX, tagP, tagXR, tagPR
        vars, nX, pars, nP, ds, orbitgrid
        ExplicitAE, HiddenAE, index, npol
        hopfX, hopfP, cputime, objval, objcount, shape
        time, solutionMatrix, eigenValues
        report, keepX = [], keepP = [], trinity, ssBranch = {}, pBranch = {}, hBranch = {}
        notSS = 0, penalty = 1 
    end
    
    properties (Access = public)
        ssopt = optimoptions(@fsolve,'Algorithm','levenberg-marquardt','Display','iter',...
                             'FunctionTolerance',1e-10,'StepTolerance',1e-10) 
                         
        dynopt = struct('abstol', 1e-10, 'reltol', 1e-8, 'solver', 'dasslc')
        optopt = struct('algorithm', 'simplex', 'tolfun', 1e-10, 'tolx', 1e-10, 'display',...
                        'iter', 'maxiter', 1e3, 'maxfunevals', 1e5, 'objectivelimitdet', 1e-9, 'nger', 100, 'npas', 100,...
                        'objectivelimitstoc', 0.0099,'objectiveshape','absolute')
                    
        conopt = struct('kmax',1000,'ds',0.01,'ncp',80,'algorithm','levenberg-marquardt','fixedDs',true);
    end        
                 
    methods (Access = public)
        %% Constructor
        function obj = Continum(model, vars, pars, massMatrix, vIndex)
            if ~exist('labudde.m','file') && isfolder('Aux')
                addpath('Aux')
            end
            if ~exist('cstr_simple.m','file') && isfolder('Models')
                addpath('Models')
            end
            obj.model = model;
            obj.ssmodel = @(x,p) model([],x,p);
            obj.vars = vars;
            obj.pars = pars;
            obj.nX = length(vars);
            obj.nP = length(pars);
            group = '[^\|\s!@#\$%&\*\(\)\-_\+={}\[\]\^~\\/\.\,;:\?ºª<>]';
            for i = obj.nX:-1:1
                obj.valueX(i,1) = vars(i).value;
                obj.lowerX(i,1) = vars(i).lower;
                obj.upperX(i,1) = vars(i).upper;
                obj.tagX{i,1} = vars(i).tag;
                obj.tagXR{i,1} = vars(i).tag(regexp(vars(i).tag, group));
            end
            if length(unique(obj.tagX)) < obj.nX
                error('Check for duplicated variable tags.')
            end
            for i = obj.nP:-1:1
                obj.valueP(i,1) = pars(i).value;
                obj.lowerP(i,1) = pars(i).lower;
                obj.upperP(i,1) = pars(i).upper;
                obj.fixedP(i,1) = pars(i).fixed;
                obj.tagP{i,1} = pars(i).tag;
                obj.tagPR{i,1} = pars(i).tag(regexp(pars(i).tag, group));
            end
            if length(unique(obj.tagP)) < obj.nP
                error('Check for duplicated parameter tags.')
            end
            obj.massMatrix = eye(obj.nX);
            if nargin > 3 && ~isempty(massMatrix)
                obj.massMatrix = massMatrix;
            end
            obj.vIndexVar = [];
            if nargin > 4 && ~isempty(vIndex)
                obj.vIndexVar = vIndex;
            end
            if exist('dasslc','file') ~= 3
                obj.dynopt.solver = 'ode15i';
            end
            obj.initialSetup;
            [err, message] = obj.check;
            if (err)
                error(message);
            end
            if obj.npol < 2
                warning('This is a first order system. No complex dynamic behavior can be observed.');
            end
            obj.reportBack('constructor');
        end
        
        %% Edition methods
%         function [p,par] = RetrieveOrbits(obj)
%            if ~isempty(obj.orbits)
%                p = obj.orbits;
%                par = obj.pBranch{2}(:,1);
%            end
%         end
        
        % Fix a parameter
        function FixParameter(obj, par)
            if strcmp(par.getType,'var')
                error('Can oly fix parameters.')
            end
            n = length(par);
            for i = 1:n
                obj.fixedP(strcmp(obj.tagP,par(i).tag)) = true;
                par(i).fixed = true;
                par(i).Show;
            end
            obj.reportBack('setfunc');
        end
        
        % Free a parameter
        function FreeParameter(obj, par)
            if strcmp(par.getType,'var')
                error('Can oly fix parameters.')
            end
            n = length(par);
            for i = 1:n
                obj.fixedP(strcmp(obj.tagP,par(i).tag)) = false;
                par(i).fixed = false;
                par(i).Show;
            end
            obj.reportBack('setfunc');
        end
        
        % Set lower bound
        function SetLower(obj, par, val)
            n = length(par);
            m = length(val);
            if n ~= m
                error('Different input lengths.')
            end
            for i = 1:n
                if strcmp(par(i).getType,'var')
                    obj.lowerX(strcmp(obj.tagX,par(i).tag)) = val(i);
                else
                    obj.lowerP(strcmp(obj.tagP,par(i).tag)) = val(i);
                end
                par(i).lower = val(i);
                par(i).Show;
            end
            obj.reportBack('setfunc');
        end
        
        % Set upper bound
        function SetUpper(obj, par, val)
            n = length(par);
            m = length(val);
            if n ~= m
                error('Different input lengths.')
            end
            for i = 1:n
                if strcmp(par(i).getType,'var')
                    obj.upperX(strcmp(obj.tagX,par(i).tag)) = val(i);
                else
                    obj.upperP(strcmp(obj.tagP,par(i).tag)) = val(i);
                end
                if par(i).value > val
                    par(i).value = val;
                end
                par(i).upper = val(i);
                par(i).Show;
            end
            obj.reportBack('setfunc');
        end
        
        % Set value
        function SetValue(obj, par, val)
            n = length(par);
            m = length(val);
            if n ~= m
                error('Different input lengths.')
            end
            for i = 1:n
                if strcmp(par(i).getType,'var')
                    obj.valueX(strcmp(obj.tagX,par(i).tag)) = val(i);
                else
                    obj.valueP(strcmp(obj.tagP,par(i).tag)) = val(i);
                end
                if par(i).value < val
                    par(i).value = val;
                end
                par(i).value = val(i);
                par(i).Show;
            end
            obj.reportBack('setfunc');
        end
        
        % Keep bifurcation point
        function KeepBP(obj)
            if isempty(obj.hopfP) 
                error('No active search')
            end
            if isempty(obj.keepX) || ~any(sum(obj.hopfP' == obj.keepP,2) == obj.nP)
                obj.keepX(end+1,:) = obj.hopfX'; 
                obj.keepP(end+1,:) = obj.hopfP'; 
            else
                warning('The active point is already listed.')
            end
            obj.reportBack('keep')
        end
        
        % Get bifurcation point
        function [x,p] = GetBP(obj, n)
            if nargin < 2
                error('Must specify a point to get.')
            elseif ~isfield(obj.report,'bplist')
                error('Must have at least one kept point.')
            elseif n > max(obj.report.bplist)
                warning('Specified point does not exist. Returning the first one.')
                n = 1;
            end
            x = obj.keepX(n,:);
            p = obj.keepP(n,:);
        end
        
        function ShiftToBP(obj, n)
            if nargin < 2
                error('Must specify a point to get.')
            elseif ~isfield(obj.report,'bplist')
                error('Must have at least one kept point.')
            elseif n > max(obj.report.bplist)
                warning('Specified point does not exist. Shifting to the first one.')
                n = 1;
            end
            obj.valueP = obj.keepP(n,:)';
            obj.valueX = obj.keepX(n,:)';
            for i = 1:obj.nX
                obj.vars(i).value = obj.valueX(i);
            end
            for i = 1:obj.nP
                obj.pars(i).value = obj.valueP(i);
            end
            obj.reportBack('setfunc')
        end
        
        %% Solver methods
        % Steady State Solver
        function [x,f,exf,out,J] = SteadyState(obj, par)
            if nargin < 2
                par = obj.valueP;
            end
            [x,f,exf,out,J] = fsolve(obj.ssmodel,obj.valueX,obj.ssopt,par);
                if exf <= 0
                    [x,f,exf,out,J] = hfsolve(obj.ssmodel,x,obj.ssopt,5,false,par);
                    if exf <= 0
                        [x,f,exf,out,J] = hfsolve(obj.ssmodel,x,obj.ssopt,25,false,par);
                    end
                end
            
            if exf > 0
                obj.valueX = x;
            elseif obj.notSS > 100
                warning('Could not find steady-state several times.')
                obj.notSS = 0;
            else
                obj.notSS = obj.notSS + 1;
            end
        end
        
        % Dynamic Solver
        function [t,y] = Simulate(obj, it, y0, yp0, par, vIndex)
            
            if nargin < 2 || isempty(it)
                error('Need integration interval.')
            end
            if nargin < 3 || isempty(y0)
                y0 = obj.valueX;
            end
            if nargin < 5 || isempty(par)
                par = obj.valueP;
            end
            if nargin < 4 || isempty(yp0)
                yp0 = obj.model(0, y0, par);
            end
            if nargin < 6
                vIndex = obj.vIndexVar;
            end
                
            opt = odeset('reltol',obj.dynopt.reltol,'abstol',obj.dynopt.abstol);
            B = obj.massMatrix;
            modelEq = obj.model;
            
            if strcmpi(obj.dynopt.solver,'ode15s') || strcmpi(obj.dynopt.solver,'ode23t')
                if rank(B) < obj.nX
                    opt.Mass = B;
                end
                [t,y] = feval(obj.dynopt.solver,obj.model, it, y0, opt, par);
            elseif strcmpi(obj.dynopt.solver,'ode15i')
                [t,y] = ode15i(@resmodel, it, y0, yp0, opt);
            elseif strcmpi(obj.dynopt.solver,'dasslc')
                [t,y] = dasslc(@resmodel, it, y0, yp0, [], obj.dynopt.reltol, obj.dynopt.abstol, vIndex);
            else
                error('Unsupported integrator.')
            end
            
            obj.time = t;
            obj.solutionMatrix = y;
            
                function [res,ires] = resmodel(~, y, yp)
                    res = B*yp - modelEq([], y, par);
                    ires = 0;
                end

        end
        
        % Hopf Search Solver
        function [x,p] = HopfSearch(obj, alg)

            function stop = outfun(~, op, ~, ~)
                if op.fval < objmin
                    stop = true;
                else
                    stop = false;
                end
            end
            
            optfmin = optimset('Display',obj.optopt.display,'TolFun',obj.optopt.tolfun,'TolX',obj.optopt.tolx,...
                                'MaxFunEvals',obj.optopt.maxfunevals,'MaxIter',obj.optopt.maxiter/5, 'Outputfcn', @outfun);
            
            objmin = obj.optopt.objectivelimitdet;                
            tic
            obj.objcount = 0;
            obj.penalty = 1;
            
            if nargin < 2
                alg = obj.optopt.algorithm;
            end
            
            swarmdisp = ~strcmpi(obj.optopt.display,'none');
            obj.valueP(obj.valueP > obj.upperP) = obj.upperP(obj.valueP > obj.upperP);
            obj.valueP(obj.valueP < obj.lowerP) = obj.lowerP(obj.valueP < obj.lowerP);
            switch lower(alg)
                
                case 'swarm'
                    [p, obj.objval, obj.objcount] = argeSwarm(@objective, obj.valueP, obj.lowerP, obj.upperP,...
                                                                      obj.optopt.nger, obj.optopt.npas, obj.optopt.objectivelimitstoc, ...
                                                                      obj.optopt.maxfunevals, swarmdisp, obj);
                case 'simplex'
                    p0 = obj.valueP;
                    for i = 1:5
                        [p, obj.objval,~,out] = fminsearch(@objective,p0,optfmin,obj);
                        obj.objcount = out.funcCount + obj.objcount;
                        p0 = p;
                        if obj.objval < obj.optopt.objectivelimitdet
                            break
                        end
                    end
                case 'hybrid'
                    
                    [p0, obj.objval, obj.objcount] = argeSwarm(@objective, obj.valueP, obj.lowerP, obj.upperP,...
                                                                      obj.optopt.nger, obj.optopt.npas, obj.optopt.objectivelimitstoc, ...
                                                                      obj.optopt.maxfunevals, swarmdisp, obj);
                    
                    for i = 1:5
                        [p, obj.objval,~,out] = fminsearch(@objective,p0,optfmin,obj);
                        obj.objcount = out.funcCount + obj.objcount;
                        p0 = p;
                        if obj.objval < obj.optopt.objectivelimitdet
                            break
                        end
                    end

                case 'levenberg-marquardt'
                    optlvm = optimoptions(@fsolve,'Algorithm',obj.optopt.algorithm,'Display',obj.optopt.display,...
                                      'TolFun',obj.optopt.tolfun,'TolX',obj.optopt.tolx,'MaxFunEvals',obj.optopt.maxfunevals,...
                                      'MaxIter',obj.optopt.maxiter);
                    [p, obj.objval,~,out] = fsolve(@objective,obj.valueP,optlvm);
                    obj.objcount = out.funcCount;
                otherwise
                    warning('Solver not recognized. Using "simplex" instead!')
                    p0 = obj.valueP;
                    for i = 1:5
                        [p, obj.objval,~,out] = fminsearch(@objective,p0,optfmin,obj);
                        obj.objcount = out.funcCount + obj.objcount;
                        p0 = p;
                        if obj.objval < obj.optopt.objectivelimitdet
                            break
                        end
                    end
            end
            obj.hopfP = obj.valueP;
            obj.hopfP(obj.fixedP == false) = p(obj.fixedP == false); 
            for i = 1:obj.nP
                obj.pars(i).value = obj.hopfP(i);
            end
            obj.valueP = obj.hopfP;
            [obj.hopfX,~,~,~,J] = obj.SteadyState;
            for i = 1:obj.nX
                obj.vars(i).value = obj.hopfX(i);
            end
            obj.cputime = toc;
            ev = eig(J,obj.massMatrix);
            obj.eigenValues = ev(~isinf(ev));
            obj.reportBack('hopfsearch')
            x = obj.hopfX; 
            p = obj.hopfP;
        end
        
        function ContinueHopf2Par(obj,par1,par2,dp)
            if nargin < 3
                error('At least two parameters must be given.')
            end
            if nargin < 4
                dp = [0.1, 0.1]';
            end
            if length(dp) < 2
                dp(2) = dp(1);
            end
            
            % Save initial state
            initialState = obj.fixedP;
            initialDisplay = obj.optopt.display;
            intialValueP = obj.valueP;
            intialValueX = obj.valueX;
            
            % Do some setup
            obj.optopt.display = 'none';
            dir = 1; n = 10000; h = 1; i = 1;
            xp = zeros(n+1,obj.nP);% f = zeros(n+1,1); 
            xp(1,:) = obj.valueP';
            pos1 = strcmp(par1.tag,obj.tagP);
            pos2 = strcmp(par2.tag,obj.tagP);
            obj.fixedP = ~strcmp(par1.tag,obj.tagP);
            
            % Calculate
            for k = 1:size(obj.keepP,1)
                obj.ShiftToBP(k)
                ok = true;
                changeDir = false;
                while ok
                    while i < n
                        if  any(obj.valueP <= obj.lowerP) || any(obj.valueP >= obj.upperP)
                            obj.valueP = lastP;
                            h = 1;
                            if changeDir
                                ok = false;
                            end
                            changeDir = true;
                            break
                        end
                        
                        [~, xp(i,:)] = obj.HopfSearch('simplex');
                        
                        if obj.objval < 1e-7
                            fprintf('iter: %4d\t fval: %2.4e\t %s: %2.4e\t %s: %2.4e\n',i,obj.objval,par1.tag,xp(i,pos1),par2.tag,xp(i,pos2))
                            lastX = obj.valueX;
                            lastP = obj.valueP;
                            obj.valueP(pos2) = obj.valueP(pos2) + dir*dp(2);                            
                            i = i+1;
                        else
                            obj.valueX = lastX;
                            obj.valueP = lastP;
                            obj.valueP(pos2) = obj.valueP(pos2) + dir*dp(2)*2^-h;
                            h = h+1;
                            if h > 10
                                h = 1;
                                if changeDir
                                    ok = false;
                                end
                                changeDir = true;
                                break
                            end
                        end
                    end
                    if changeDir
                        obj.ShiftToBP(k)
                        dir = -dir;
                        xp(i,:) = nan;
                        i = i+1;
                    end
                end
            end
            xp(i:end,:) = [];
            obj.fixedP = initialState;
            obj.valueP = intialValueP;
            obj.valueX = intialValueX;
            obj.optopt.display = initialDisplay;
            
            obj.hBranch{end+1} = par1.tag;
            obj.hBranch{end+1} = par2.tag;
            obj.hBranch{end+1} = xp;
        end
        
        function PeriodicContinuation(obj,par,vIndex,vars)
            
            if nargin < 2
                error('One parameter must be given')
            end
            if nargin < 3 || isempty(vIndex)
                vIndex = obj.vIndexVar;
            end
            if nargin < 4
                vec = 1:obj.nX;
            end
            if nargin == 4
                nv = length(vars);
                vec = zeros(nv,1);
                for j = 1:nv
                    posj = find(strcmp(vars(j).tag,obj.tagX));
                    if ~isempty(posj)
                        vec(j) = posj;
                    else
                        error('Problem does not contain variable "%s".', vars(j).tag)
                    end
                end
            end
            
            if ~ismember(obj.hopfP',obj.keepP,'rows')
                obj.ShiftToBP(1)
            end
            
            % The problem function
            function F = theProblem(x)
                
                % Alocate the "variables"
                Y = reshape(x(1:n*m),n,m);
                T = x(n*m+1);
                if ~initialOrbitSearch
                    newpar(pos) = x(n*m+2);
                end
                
                % Discretization
                for i = m:-1:1
                    Fk(:,i) = cModel([], Y(:,i), newpar);
                end
                R = (B*Y)*Lt - T*Fk;
                
                % The equation's residuals
                R = R(n+1:end)';
                
                % The boundary condition
                P = Y(:,1) - Y(:,end);
                
                % The inflated integral
                I1 = theta(1)*trapz(tgrid,sum((Y(vec,:)-Y0(vec,:)).*Ydot(vec,:))) + theta(2)*(T-T0)*Tdot + theta(3)*(newpar(pos) - par0)*pardot - deltaS; 
                                
                % The phase condition
                if ~initialOrbitSearch
                    I2 = trapz(tgrid,sum(Y(vec,:).*dY0(vec,:)));       
                    F = [R; I1; I2; P];
                else
                    F = [R; I1; P];
                end

            end
            
            % The initial parameters
            B = obj.massMatrix;
            n = obj.nX;
            m = obj.conopt.ncp;
            cModel = obj.model;
            newpar = obj.valueP;
            deltaS = obj.conopt.ds;
            ds0 = deltaS;
            pos = strcmp(par.tag,obj.tagP);
            theta = [1 1 1];
            
            % The collocation stuff
            [L, ~, tgrid] = collocation(m-2, '2*x-1', [0 1]);
            Lt = L';
            obj.orbitgrid = tgrid;
            
            % Solve the problem at current HBP
            [xi,~,~,~,J] = obj.SteadyState(newpar);
            
            % Find the initial period
            evals = eig(J,B);
            evals(isinf(evals)) = [];
            epair = evals(abs(real(evals)) < 1e-6 & abs(imag(evals)) > 1e-7);
            T0 = 2*pi/unique(abs(imag(epair)));
            
            % Solver setup
            opt = obj.ssopt;
            opt.Display = 'iter';
            opt.Algorithm = obj.conopt.algorithm;
            opt.MaxFunEvals = 2e5;
            opt.MaxIter = 2e4;
            
            % Find the initial orbit
            initialOrbitSearch = true;
            
            % Problem function parameters
%             Ydot = zeros(n,m);
%             for j = 1:m
%                 Ydot(:,j) = expm((J*T0)*tgrid(j))*xi; % change to a general analytic solution. Now -> ODE only
%             end
            Y0 = xi.*ones(n,m);
            Y0a = Y0;
            Tdot = 0;
            pardot = 0;
            par0 = newpar(pos);
            
            function [res,ires] = linprob_fi(~,y,yp)
                res = B*yp - J*y;
                ires = 0;
            end
            
            function dy = linprob(~,y)
                dy = J*y;
            end
            
            if ~isempty(vIndex)
                xpi = -linprob_fi(0,xi,zeros(obj.nX,1));
                [~,Ydot2] = dasslc(@linprob_fi,tgrid*T0, xi, xpi, [], obj.dynopt.reltol, obj.dynopt.abstol, vIndex);
            else
                [~,Ydot2] = ode15s(@linprob,tgrid*T0, xi, odeset('abstol',obj.dynopt.abstol,'reltol',obj.dynopt.reltol,'Mass',B));
            end
            Ydot = Ydot2';
            
            % The initial guess
            x0 = xi.*ones(n,m);
            x0 = reshape(x0,n*m,1);
            x0(n*m+1) = T0;
 
            % Solver call
            fprintf('Searching for initial orbit...\n\n')
            [x,f] = fsolve(@theProblem, x0, opt);
           
            if norm(f) > 1e-4
                error('No initial orbit detected.')
            end

            % The results
            po = reshape(x(1:n*m),n,m);
            
            % Begin continuation
            initialOrbitSearch = false;
            
            % Insert the continuation parameter in the variable's vector
            x(n*m+2) = newpar(pos);
            
            % Continuation parameters
            kend = obj.conopt.kmax;
            parspace = zeros(kend,1);
            maxsol = zeros(kend,obj.nX);
            parspace(1) = newpar(pos);
            maxsol(1,:) = max(po,[],2)';
            opt.Display = 'none';
            orbits = zeros(n,m,kend);
            orbits(:,:,1) = po;
            periods = zeros(kend,1);
            periods(1) = x(n*m+1);
            floquet = zeros(obj.npol,kend);
            [~,~,~,~,Jss] = obj.SteadyState();
            [V,U,i1] = pencan(obj.massMatrix,Jss);
            KCF = U*Jss*V;
            KCF_J = KCF(1:i1,1:i1);
            floquet(:,1) = eig(expm(KCF_J*x(n*m+1)));
            
            iterate = true;
            k = 2;
            cont = 1;
            sf = sprintf('\nRuning periodic continuation...\n\n');
            fprintf([sf repmat(' ',1,50)])
            
            while iterate
                Y0 = po;
                dY0 = (Y0(:,2:end) - Y0(:,1:end-1))./(tgrid(2:end)' - tgrid(1:end-1)');
                dY0(:,m) = dY0(:,m-1);
                Tdot = 1/deltaS*(x(n*m+1) - T0);
                pardot = 1/deltaS*(x(n*m+2) - par0);
                T0 = x(n*m+1);
                par0 = x(n*m+2);
                Ydot = 1/deltaS*(Y0 - Y0a);
                Y0a = Y0;
%                 theta = 1./[norm(Ydot)^2 Tdot^2 pardot^2+sqrt(eps)];
                [x,f] = fsolve(@theProblem, x, opt);
                po = reshape(x(1:n*m),n,m);
                if norm(f) < 1e-8
                    parspace(k) = x(n*m+2);
                    maxsol(k,:) = max(po,[],2)';
                    orbits(:,:,k) = po;
                    periods(k) = x(n*m+1);
                    [~,~,~,~,Jss] = obj.SteadyState();
                    [V,U,i1] = pencan(obj.massMatrix,Jss);
                    KCF = U*Jss*V;
                    KCF_J = KCF(1:i1,1:i1);
                    floquet(:,k) = eig(expm(KCF_J*x(n*m+1)));
                    sf = sprintf('Iter: %4d\t par: %2.6e\t norm: %2.6e\n',k,parspace(k),norm(f));
                    fprintf([repmat('\b',1,50) sf])
                    k = k + 1;
                    cont = 1;
                    if ~obj.conopt.fixedDs
                        deltaS = min(deltaS+ds0, 8*ds0);
                    end
                    if k > kend
                        iterate = false;
                    end
                else
                    if ~obj.conopt.fixedDs
                        deltaS = max(deltaS/2, ds0/10);
                    end
                    cont = cont + 1;
                    if cont > 10
                        warning('Could not converge in continuation procedure.')
                        sf = sprintf('Iter: %4d\t par: %2.6e\t norm: %2.6e\n',k-1,parspace(k-1),norm(f));
                        fprintf([repmat('\b',1,50) sf])
                        iterate = false;
                    end
                end
                
            end
            fprintf('\nFinished!')
            if k < kend
                parspace(k:end) = [];
                maxsol(k:end,:) = [];
                orbits(:,:,k:end) = [];
                periods(k:kend) = [];
                floquet(:,k:kend) = [];
            end
            if isempty(obj.pBranch)
                obj.pBranch{1} = par.tag;
                obj.pBranch{2} = [parspace, maxsol];
                obj.pBranch{3} = orbits;
                obj.pBranch{4} = periods;
                obj.pBranch{5} = floquet;
            else
                obj.pBranch{end+1} = par.tag;
                obj.pBranch{end+1} = [parspace, maxsol];
                obj.pBranch{end+1} = orbits;
                obj.pBranch{end+1} = periods;
                obj.pBranch{end+1} = floquet;
            end

        end
        
        function SSContinuation(obj,par,n,it)
            
            if nargin < 2
                error('You need to specify the parameter.')
            end
            if length(par)>1
                error('Only one parameter is allowed.')
            end
            if nargin < 3 || isempty(n)
                n = 1000;
            end
            if nargin < 4 
                it = [par(1).lower, par(1).upper];
            end
            
            B = obj.massMatrix;
            parspace = linspace(it(1),it(2),n)';
            pos = strcmp(par.tag,obj.tagP);
            x = zeros(n,obj.nX);
            lambda = zeros(n,obj.nX);
            newpar = obj.valueP;
            sf = sprintf('\nRuning steady-state continuation...\n\n');
            fprintf([sf repmat(' ',1,31)])
            for i = 1:n
                newpar(pos) = parspace(i);
                [x(i,:),f,~,~,J] = obj.SteadyState(newpar);
                lambda(i,:) = eig(J,B);
                sf = sprintf('Iter: %4d\t norm: %2.6e\n',i,norm(f));
                fprintf([repmat('\b',1,31) sf])
            end
            fprintf('\nFinished!')
            r = real(lambda);

            stableC = nan(n,obj.nX);
            unstableC = stableC;
            for i = 1:n
                if any(r(i,~isinf(r(i,:)))>0)
                    unstableC(i,:) = x(i,:);
                else
                    stableC(i,:) = x(i,:);
                end
            end

            if isempty(obj.ssBranch)
                obj.ssBranch{1} = par.tag;
                obj.ssBranch{2} = [parspace, stableC, unstableC]; 
                obj.ssBranch{3} = lambda;
            else
                obj.ssBranch{end+1} = par.tag;
                obj.ssBranch{end+1} = [parspace, stableC, unstableC];
                obj.ssBranch{end+1} = lambda;
            end
        end

        %% Report method
        function ShowReport(obj)
            fprintf('\n\n==================== Reports ====================\n\n')
            fprintf('* Initial Setup\n\n')
            fprintf('    Model ........................................... %s\n',[func2str(obj.model) '.m']);
            fprintf('    Variables ....................................... ')
            for i = 1:obj.nX
                fprintf('%s = %d, ', obj.tagXR{i}, obj.report.x0(i));
            end
            fprintf('\b\b\n') 
            fprintf('    Search Parameters ............................... ')
            for i = 1:obj.nP
                if ~obj.fixedP(i)
                    fprintf('%s = %d, ', obj.tagPR{i}, obj.report.par0(i));
                end
            end
            fprintf('\b\b\n')
            if any(obj.fixedP)
                fprintf('    Fixed Parameters ................................ ')
                for i = 1:obj.nP
                    if obj.fixedP(i)
                        fprintf('%s = %d, ', obj.tagPR{i}, obj.report.par0(i));
                    end
                end
                fprintf('\b\b\n\n')
            else
                fprintf('\n')
            end
            fprintf('* Model Analysis\n\n')
            fprintf('    Differential Index .............................. %d\n', obj.report.DAE_Index);
            fprintf('    Explicit Algebraic Eqs .......................... %d\n', obj.report.ExplicitAE);
            fprintf('    Hidden Algebraic Eqs ............................ %d\n', obj.report.HiddenAE);
            fprintf('    Dynamic Degrees of Freedom ...................... %d\n\n', obj.report.DDoF);

            if isfield(obj.report,'x')
                fprintf('* Active Search\n\n')
                fprintf('    Hopf found? ..................................... ')
                if obj.report.objval <= obj.optopt.objectivelimitdet
                    fprintf('Yes.\n\n')
                elseif obj.report.objval > obj.optopt.objectivelimitdet && obj.report.objval <= obj.optopt.objectivelimitdet*100
                    fprintf('Yes. But more refinment may be needed.\n\n')
                else
                    fprintf('No.\n\n')
                end
                fprintf('    x-values ........................................ ')
                for i = 1:obj.nX
                    fprintf('%s = %d, ', obj.tagXR{i}, obj.report.x(i));
                end
                fprintf('\b\b\n');
                fprintf('    par-values ...................................... ')
                for i = 1:obj.nP
                    if ~obj.fixedP(i)
                        fprintf('%s = %d, ', obj.tagPR{i}, obj.report.par(i));
                    end
                end
                fprintf('\b\b\n');
                fprintf('    Eigenvalues ..................................... (')
                for i = 1:obj.npol
                    fprintf('%.4e %+.4ei, ', real(obj.report.eigenValues(i)), imag(obj.report.eigenValues(i)))
                end
                fprintf('\b\b)\n\n');
                fprintf('    CPU Time ........................................ %d s\n', obj.report.cputime);
                fprintf('    Obj Function Value .............................. %d\n', obj.report.objval);
                fprintf('    Function Count .................................. %d\n\n', obj.report.objcount);
            end
            
            if isfield(obj.report,'keepX')
                fprintf('* Kept Bifurcation Points\n\n')
                for i = 1:length(obj.report.bplist)
                    fprintf('    %d:    x = (',i)
                    fprintf('%d, ', obj.report.keepX(i,:))
                    fprintf('\b\b)\n')
                    fprintf('    \t  p = (')
                    fprintf('%d, ', obj.report.keepP(i,:))
                    fprintf('\b\b)\n\n')
                end
            end
            fprintf('\n================= End of Report =================\n\n')

        end
        
        %% Plot methods
        function PlotTimeSeries(obj,var)
            if nargin < 2
                error('Must specify a variable')
            end
            
            if any(~strcmp(var.getType,'var'))
                error('Only variables can be plotted.')
            end
            
            for i = length(var):-1:1
                pos(1,i) = find(strcmp(obj.tagX,var(i).tag));
            end
            
            figured;
            plot(obj.time, obj.solutionMatrix(:,pos),'linewidth',1.5)
            xlabel('Time')
            if length(var) > 1
                legend(var.tag,'Location','Best')
                ylabel('Variables')
            else
                ylabel(var.tag)
            end
        end
        
        function PlotPhasePlane(obj,var1,var2)
            
            if nargin < 3
                error('Must specify a pair of variables.')
            end
            
            if length(var1) > 1 || length(var2) > 1
                error('Both inputs must be scalars.')
            end
            
           
            figured;
            plot(obj.solutionMatrix(:,strcmp(obj.tagX,var1.tag)),obj.solutionMatrix(:,strcmp(obj.tagX,var2.tag)),'linewidth',1.5)
            hold on
            plot(obj.solutionMatrix(1,strcmp(obj.tagX,var1.tag)),obj.solutionMatrix(1,strcmp(obj.tagX,var2.tag)),...
                's','MarkerFaceColor',[0.85,0.325,0.098],'MarkerSize',12,'MarkerEdgeColor',[0.85,0.325,0.098])
            plot(obj.solutionMatrix(end,strcmp(obj.tagX,var1.tag)),obj.solutionMatrix(end,strcmp(obj.tagX,var2.tag)),...
                'o','MarkerFaceColor',[0.929,0.694,0.125],'MarkerSize',10,'MarkerEdgeColor',[0.929,0.694,0.125])
            xlabel(var1.tag)
            ylabel(var2.tag)
            legend('Trajectory','Start','End','Location','Best')
            
        end
        
        function PlotEigenMap(obj, par, n, it)
            if nargin < 2
                error('You need to specify a parameter.')
            end
            if nargin < 3 || isempty(n)
                n = 100;
            end
            if nargin < 4 
                it = [par.lower, par.upper];
            end
            
            parspace = linspace(it(1),it(2),n)';
            pos = strcmp(par.tag,obj.tagP);
            newpar = obj.valueP;
            for i = n:-1:1
                newpar(pos) = parspace(i);
                [obj.valueX,~,~,~,Jss] = obj.SteadyState(newpar);
                eigs = eig(Jss,obj.massMatrix)';
                lambda(i,:) = eigs(~isinf(eigs));
            end
            realL = real(lambda);
            complexL = imag(lambda);
            
            obj.trinity = [parspace, realL, complexL];
            
            newFig = figured;
            plot(realL,complexL,'LineWidth',1,'LineStyle','-','Marker','.','MarkerSize',10);
            ax = gca;
            line([0,0],ax.YLim,'color','k','linestyle',':')
            line(ax.XLim,[0,0],'color','k','linestyle',':')
            set(gca,'YGrid','on','Box','on','FontSize',18)
            ylabel('Im(\lambda)','FontWeight','bold','FontName','Palatino','fontsize',18)
            xlabel('Re(\lambda)','FontWeight','bold','FontName','Palatino','fontsize',18)
    
            
            cursorMode = datacursormode(newFig);
            hPlot = get(gca,'Children');

            nLines = length(hPlot);
            hDatatip{nLines-2,1} = 0;
            XData{nLines-2,1} = 0;
            YData{nLines-2,1} = 0;
            set(cursorMode,'UpdateFcn',{@eigenCursor,[parspace, lambda]},'DisplayStyle','datatip')
            k = 1;
            rec = 0;
            for i = 1:nLines
                if length(hPlot(i).XData) > 2
                    XData{k} = hPlot(i).XData;
                    YData{k} = hPlot(i).YData;
                    hDatatip{k} = cursorMode.createDatatip(hPlot(i));
                    set(hDatatip{k}, 'Marker','o', 'MarkerSize',10, 'MarkerFaceColor',hPlot(i).Color,...
                             'MarkerEdgeColor','w', 'OrientationMode','auto',...
                             'Position', [hPlot(i).XData(1) hPlot(i).YData(1) 0])
                    updateDataCursors(cursorMode)
                    drawnow;
                    k = k+1;
                end
            end
            set(gca,'FontSize',16,'Position',[0.13 0.11 0.775 0.815])
            rec = 1;
            ind = 1;
        %   
                function txt = eigenCursor(~,event_obj,data)

                    % cursorMode = datacursormode(gcf);
                    pos = get(event_obj,'Position');
                    I = get(event_obj, 'DataIndex');
                    txt = {['Re : ',sprintf('% 2.4g',pos(1))],...
                           ['Im : ',sprintf('% 2.4g',pos(2))],...
                           ['Par: ',sprintf('% 2.8g',data(I,1))]};

                    if rec == 1
                        rec = 0;
                        ind = I;
                        for j = 1:nLines - 2
                            hDatatip{j}.Cursor.Position = [XData{j}(ind) YData{j}(ind) 0];
                            drawnow;
                        end
                        rec = 1;
                    end

                end
        end
        
        function PlotObjectiveFunction(obj, par, n, it)
            if nargin < 2
                error('You need to specify the parameter.')
            end
            if length(par)>2
                error('Maximum of two parameters is allowed.')
            end
            if nargin < 3 || isempty(n)
                n = [100 100];
            end
            if length(n) < 2
                n(2) = n(1);
            end
            if nargin < 4 
                it = [par(1).lower, par(1).upper; par(2).lower, par(2).upper];
            end
            
            obj.penalty = 0; %Becareful
            switch length(par)
                case 1
                    parspace = linspace(it(1,1),it(1,2),n(1))';
                    pos = strcmp(par.tag,obj.tagP);
                    newpar = obj.valueP;
                    for i = n(1):-1:1
                        newpar(pos) = parspace(i);
                        f(i) = objective(newpar, obj);
                    end
                    
                    figured;
                    plot(parspace,f,'LineWidth',1.5);
                    xlabel(par.tag)
                    ylabel('Objective function')
                case 2
                    par1space = linspace(it(1,1),it(1,2),n(1))';
                    par2space = linspace(it(2,1),it(2,2),n(2))';
                    pos = logical(strcmp(par(1).tag,obj.tagP) + strcmp(par(2).tag,obj.tagP));
                    newpar = obj.valueP;
                    for i = n(1):-1:1
                        for j = n(2):-1:1
                            newpar(pos) = [par1space(i) par2space(j)];
                            f(j,i) = objective(newpar, obj);
                        end
                    end
                    [P1,P2] = meshgrid(par1space,par2space);
                    figured;
                    surfc(P1,P2,f,'EdgeColor','none');
            end
            obj.penalty = 1;
        end
        
        function PlotBranchDiagram(obj,par,var,points)
            
            if nargin < 3
                error('You need to specify a parameter and a variable.')
            end
            if nargin < 4 || isempty(points)
                points = inf;
            end
            
            txt = {};
            holdon = false; 
            
            if isempty(obj.ssBranch) || ~any(strcmp(obj.ssBranch,par.tag))
                warning('No stationary branches detected for the parameter "%s". \n Try calling "SSContinuation" method first.', par.tag)
            else
                posX = strcmp(var.tag,obj.tagX);
                posP = strcmp(par.tag,obj.tagP);
                posX2 = find(posX);
                figured;
                hold on
                for i = 1:3:length(obj.ssBranch)
                    if strcmp(obj.ssBranch{i},par.tag)
                        parspace = obj.ssBranch{i+1}(:,1);
                        stableC = obj.ssBranch{i+1}(:,1+posX2);
                        unstableC = obj.ssBranch{i+1}(:,1+obj.nX+posX2);
                        p1 = plot(parspace, stableC,'k-',parspace,unstableC,'k:');
                        set(p1(1:2),'LineWidth',1.5)
                        xlabel(par.tag)
                        ylabel(var.tag)
                        [~,~,~,txt] = legend(p1(1:2),{'Stable steady-state solutions', 'Unstable steady-state solutions'});
                    end
                end
                hold off
                if ~isempty(obj.keepX)
                    hold on
                    hP = obj.keepP(:,posP);
                    hX = obj.keepX(:,posX);
                    p2 = plot(hP,hX,'ks','MarkerFaceColor','k','MarkerSize',10);
                    txt{end+1} = 'Hopf bifurcation point';
                    legend([p1(1:2); p2(end)], txt)
                    hold off
                end
            end
            
            if isempty(obj.pBranch) || ~any(strcmp(obj.pBranch,par.tag))
                warning('No periodic branches detected for the parameter "%s". \n Try calling "PeriodicContinuation" method first.', par.tag)
            else
                if isempty(txt)
                    figured;
                    xlabel(par.tag)
                    ylabel(var.tag)
                else
                    hold on
                    holdon = true;
                end
                posX = strcmp(var.tag,obj.tagX);
                for i = 1:5:length(obj.pBranch)
                    if strcmp(obj.pBranch{i},par.tag)
                        parspace = obj.pBranch{i+1}(:,1);
                        maxsol = obj.pBranch{i+1}(:,logical([0; posX]));
                        last = [0; 0];
                        for j = 1:min(points,length(maxsol))
                            actual = [parspace(j); maxsol(j)];
                            if abs(last(1) - actual(1)) > abs(max(parspace) - min(parspace))/12 || ...
                                    abs(last(2) - actual(2)) > abs(max(maxsol) - min(maxsol))/12
                                if any(abs(obj.pBranch{i+4}(:,j)) > 1+1e-6)
                                    mark = 'w';
                                else
                                    mark = 'k';
                                end
                                p3 = plot(parspace(j), maxsol(j), 'ko', 'MarkerFaceColor', mark);
                                last = actual;
                            end
                        end
                    end
                end
            end
            
            txt{end+1} = 'Periodic solution (max value)';
            if holdon    
                legend([p1(1:2);p2;p3],txt,'location','best')
                hold off
            else
                legend(txt,'location','best')
            end
        end
        
        function PlotOrbitSurface(obj, par, var1, var2, step, cmap)
            
            if nargin < 4
                error('You need to specify a parameter and two variables.')
            end
            if nargin < 5 || isempty(step)
                step = 1;
            end
            if nargin < 6 || isempty(cmap)
                cmap = @parula;
            end
            
            if isempty(obj.pBranch) || ~any(strcmp(obj.pBranch,par.tag))
                warning('No orbits detected for the parameter "%s". \n Try calling "PeriodicContinuation" method first.', par.tag)
            else
                pos1 = strcmp(var1.tag,obj.tagX);
                pos2 = strcmp(var2.tag,obj.tagX);
                lo = ones(1,length(obj.orbitgrid));
                figured;
                xlabel(var1.tag)
                ylabel(var2.tag)
                zlabel(par.tag)
                hold on
                for i = 1:5:length(obj.pBranch)
                    if strcmp(obj.pBranch{i},par.tag)
                        parspace = obj.pBranch{i+1}(:,1);
                        lp = length(parspace);
                        color = cmap(lp);
                        maxbase = 0;
                        minbase = inf;
                        p = gobjects(lp,1);
                        for j = 1:step:lp
                            p(j) = plot3(obj.pBranch{i+2}(pos1,:,j),obj.pBranch{i+2}(pos2,:,j), ...
                                parspace(j)*lo,'color',color(j,:),'linewidth',1.5);
                            arclen = sum(sqrt(1+diff(obj.pBranch{i+2}(pos1,:,j)).^2));
                            if arclen > maxbase
                                maxbase = arclen;
                            end
                            if arclen < minbase
                                minbase = arclen;
                            end
                        end
                        for j=1:step:lp
                            fat = (sum(sqrt(1+diff(obj.pBranch{i+2}(pos1,:,j)).^2))-minbase)/(maxbase-minbase);
                            c= (color(end,:)*fat + (1-fat)*color(1,:));
                            p(j).Color = c;
                        end
                    end
                end
                view([45 25])
                hold off
            end
            
        end
        
        function Plot2ParDiagram(obj)
            
            if isempty(obj.hBranch)
                warning('No data detected. Try calling "Continue2ParHopf" method first.')
            else
                for i = 1:3:length(obj.hBranch)
                    figured;
                    pos1 = strcmp(obj.hBranch{i},obj.tagP);
                    pos2 = strcmp(obj.hBranch{i+1},obj.tagP);
                    plot(obj.hBranch{i+2}(:,pos1),obj.hBranch{i+2}(:,pos2),'linewidth',1.5,'color','k')
                    xlabel(obj.hBranch{i})
                    ylabel(obj.hBranch{i+1})
                end
            end
            
        end
    end
    
    methods (Access = private)
        %% Check if problem is well-posed
        function [err, message] = check(obj)
            dy = obj.model([],obj.valueX,obj.valueP);
            err = true;
            if any(isnan(dy))
                message = 'The model return NaN at initial guess.';
            elseif length(dy) ~= obj.nX
                message = sprintf('Model returns %d equations, but there are %d variables.', [length(dy), obj.nX]);
            elseif any(obj.lowerP >= obj.upperP)
                message = 'Inconsistent bounds.';
            elseif diff(size(obj.massMatrix))
                message = 'Mass matrix must be square.';
            elseif length(dy) ~= length(obj.massMatrix)
                message = 'Dimensions of mass matrix and system must agree.';
            elseif size(dy,2) > 1
                message = 'Model must return a column vector.';
            else
                err = false;
                message = 'ok';
            end
        end
        
        %% Initial setup
        function initialSetup(obj)
            I = eye(obj.nX);
            x0 = obj.valueX; %rand(obj.nX,1);
            f0 = obj.ssmodel(x0,obj.valueP);
            delta = sqrt(eps);
            for m = obj.nX:-1:1
                J_ss0(:,m) = ( obj.ssmodel(x0+delta*I(:,m),obj.valueP) - f0 )/delta;
            end
            if rcond(J_ss0) < delta
                disp0 = obj.ssopt.Display;
                valX0 = obj.valueX;
                obj.ssopt.Display = 'none';
                [~,~,~,~,J_ss0] = obj.SteadyState();
                obj.ssopt.Display = disp0;
                obj.valueX = valX0;
            end
            r = rank(obj.massMatrix);
            L0 = eig(J_ss0,obj.massMatrix);
            nPol = length(L0(~isinf(L0)));

            obj.ExplicitAE = (obj.nX - r);
            obj.HiddenAE = (r - nPol);
            [Q,M,i1] = pencan(obj.massMatrix,J_ss0);
            if i1 == obj.nX
                obj.index = 0;
            else
                KCF = M*obj.massMatrix*Q;
                N = KCF(i1+1:end,i1+1:end);
                for m = 1:obj.nX-i1
                    if norm(N^m) < sqrt(eps)%all(abs(N^m) < 1e-12)
                        obj.index = m;
                        break
                    end
                end
            end
            obj.npol = nPol;
            
            if strcmp(obj.optopt.objectiveshape,'quadratic')
                obj.shape = 1;
            else
                obj.shape = 0;
            end
        end
        
        %% Build reports
        function reportBack(obj, caller)
            switch caller
                case 'constructor'
                    obj.report.model = obj.model;
                    obj.report.x0 = obj.valueX;
                    obj.report.par0 = obj.valueP;
                    obj.report.DAE_Index =  obj.index;
                    obj.report.ExplicitAE = obj.ExplicitAE;
                    obj.report.HiddenAE = obj.HiddenAE;
                    obj.report.DDoF = obj.npol;
                case 'hopfsearch'
                    obj.report.x = obj.hopfX;
                    obj.report.par = obj.hopfP;
                    obj.report.cputime = obj.cputime;
                    obj.report.objval = obj.objval;
                    obj.report.objcount = obj.objcount;
                    obj.report.eigenValues = obj.eigenValues;
                case 'keep'
                    obj.report.keepX = obj.keepX;
                    obj.report.keepP = obj.keepP;
                    obj.report.bplist = (1:size(obj.keepX,1))';
                case 'setfunc'
                    for i = 1:obj.nX
                        obj.report.x0(i) = obj.vars(i).value;
                    end
                    for i = 1:obj.nP
                        obj.report.par0(i) = obj.pars(i).value;
                    end
                    if isfield(obj.report,'x')
                        obj.report = rmfield(obj.report,{'x'});
                    end
            end
        end
        
    end
    
        
end

%% Aux functions
function f = objective(par,obj)
    
    par(obj.fixedP==true) = obj.valueP(obj.fixedP==true);
    obj.ssopt.Display = 'none';
    [x_ss,f_ss,status,~,J_ss] = obj.SteadyState(par);
    if status <= 0
        f = 1000*pi;
        return
    end

    if norm(f_ss) > 1e-6 || any(imag(x_ss)~=0) || any(any(imag(J_ss)~=0))
        %warning('Steady-state not found!');
        f = 10000*pi;
        return
    end

    [~,P] = labudde(J_ss,obj.massMatrix);
    nPol2 = length(P) - 1;
    if obj.npol ~= nPol2
        warning('Order changed!');
        fprintf('\nStarting n = %d and Current n = %d\n',nPol,nPol2);
        f = 10000*pi/2;
        return
    end

    % Build Hurwitz matrix
    for j = nPol2-1:-1:1
        for k = 1:nPol2-1
            pos = (2*j-k+1);
            if pos > nPol2+1 || pos < 1 
                Delta(j,k) = 0;
            else
                Delta(j,k) = P(pos);
            end
        end
    end

    con = -min(0,P(end));
    bar = -min(0,min(par - obj.lowerP)) - min(0,min(obj.upperP - par));

    f = (norm(f_ss) + pi*1e6*bar + pi*1e6*con)*obj.penalty + ( det(Delta)^2 )^(1/2 + obj.shape/2);
end
