% % MAIN2.M
% %
% % USAGE: Main script for solving and analyzing sequential equilibria in 
% % our model.
% %
% % =======================================================================
% %     (c) 2013-- T.Kam (Revision for JEDC)
% %
% % DEPENDENCIES:
% %
% %     * KSMOD3 class for model and SSE solution.
% %     * GLPKMEX and GNU Linear Programming Kit (ANSI C to MATLAB MEX)
% %     * YALMIP (for BMIBNB global optimization)
% %     * SNOPT (Local nonlinear optimizer)
% %     * CompEcon Toolbox (Miranda-Fackler)
% %         
% % Email: tcy.kam__at__gmail.com
% % =======================================================================
% % $Revision: 5.0.0 $  $Date: 2013/11/08 23:38:20 $  

clear all
close all
clc

%% STARTUP RUN OPTIONS
LOAD_OLD = 'yes';               % load W_old1.mat
iter_last = 5;                  % Check number in /_mat folder

%% MACHINE SETTINGS

    % Determine number of parallel workers:
    MaxNumCores = 16;
    POOLSIZE = min( feature('numCores'), MaxNumCores );
    
    
%% PROGRAM SETTINGS

    DISPLAY_DETAIL_INNER = 0;       % Display option for while/for loops

    USE_GLPK = 1;                   % Option: GLPK or LINPROG (in-built)
    USE_PARFOR = 1;                 % Run loops in parallel
    USE_PROFILER = 0;               % Development: profile bottlenecks

    % Pack these options in object:
    Options_Sol.glpk = USE_GLPK;
    Options_Sol.pfor = USE_PARFOR;
    Options_Sol.disp = DISPLAY_DETAIL_INNER;
    
    % Plot all partition elements of state space and random elements of
    % intersect( P(a)(Q_j), Q_k) ), for all j,k = 1,...,K.
    
    plot_partsim = 'off';      % Show realizations of HnR rv's in partitions
                              % Q(k) under each mapping P(a)(Q(j))
    plot_partsiminverse = 'off';  % Actual rv's in each Q(j)
    
    % Run checks on current state:
    if USE_PROFILER == 1
        profile on -history
    end
    
    %if USE_PARFOR == 1
    %     if matlabpool('size') == 0
    %        matlabpool(POOLSIZE);
    %    end
   % end
    

    fprintf('MAIN.M: DIRTY HACK: Using iter101 from TEST_MAIN.M !\n');
    % Initial guess of hyperplane levels (payoff correspondence W_outer)
    %load sse_005					% W_old
	%load('_mat/sse2014b_init0.mat');
    
%% MAIN SECTION FOR APPROXIMATE BILINEAR FORMULATION OF SPE VALUE SETS 

if strcmp(LOAD_OLD,'no') == 1

    fprintf('MAIN.M: Starting from scratch!\n');
    
    %% STEP 0: MODEL PARAMETERS

    SIGMA = 1;                % CRRA parameter (linear = 0; log = 1)
    PHI = 0;                  % Disutility of effort parameter (linear = 0)
    M = 1;                    % "Max duration" of employment, M
    N = 2;                    % "Max duration" of unemployment, N
    WAGE = 1;                 % fixed common wage for employed
    REP_RATIO = 0.8;          % "Maximal replacement ratio"
    DELTA = 0.9615;           % common discount factor: Annual rate = 4%

    L = 192;                  % # of vertices per state partition element
    %NGRID_P = 8;             % Number of grid elements representing [0,1]
    K = 16;                   % # partition elements Q(k) of state space
    
    % NOTES: Dependant parameters set in KSMOD4.M
    % -------------------------------------------
    % MBAR = WAGE*REP_RATIO;    % exogenous upperbound on benefit payout
    % DELSTAR = (1 - DELTA);    % For averaging payoffs

    %% STEP 1: CREATE INSTANCE OF CLASS: KSMOD04
    
    fprintf('\nCreating MODEL instance from KSMOD04 ... Please wait.\n');

    model = ksmod04(SIGMA,PHI,M,N,WAGE,REP_RATIO,DELTA); 

  

    %% STEP 2: PRELIMINARIES
    
    fprintf('\nStarting Preliminary calculations ... Please wait\n');
    
    % Monte Carlo elements
    HitnRun = 0;               % Use Hit-and-Run algorithm
    Nsim = 100;                 % #uniform r.v. in each simplex intersection
    
    % Deterministic Elements
    N_Z = model.NS;                               % Private state dimension
    N_BGRID = size(model.ProfileSetB,1);              % Grid on policy b(j)
                                                  % where j = 1, ...,N_Z

    %np = N_Z + 1;                                % Payoff vector dimension
                                                  % +1 element: govt payoff   
    np = N_Z;
    
    NA = model.NA;
    
    %% Fix search subgradients. Use Leopardi (2006):      
    fprintf('1. Generate set of spherical codes H. ...\n');
    
    H = subgradient(np,L);                      % H is finite subset of the 
                                                % np-dimensional sphere: 
                                                % S^(np-1). So, for each i,
                                                % sum( H(i,:).^2, 2 ) = 1.

    fprintf('\t\t\t\t ... DONE!\n');
    
    
    %% Compute natural bounds on actions and payoffs:
    fprintf('2. Compute natural payoff bounds: V_MIN, V_MAX. ...\n');
    
    [V_MIN, V_MAX, ~] = payoff_bounds(model);
    
    % Lower/Upper bounds for GLPK later
    lb = V_MIN*ones(np, 1); % Lower bound on payoff vector
    ub = V_MAX*ones(np, 1); % Upper bound on payoff vector

    fprintf('\t\t\t\t ... DONE!\n');
    
    %% Space of private agents' action profiles: 
    fprintf('3. Preallocate (finite) action profile set, A. ...\n');
    
    %[A, NA, N, A_u, A_e] = A_ProfileSet(model);  
    
    fprintf('\t\t\t\t ... DONE!\n');
    
    %% Cartesian product space of A x B in terms of indexes:

    IndexProfileSet = ProfileSetIndex(model); % Col#1: a, and Col#2: b

    %% Markov transition law at each action profile:
    fprintf('4. Map action profile set to Markov laws, P(a)....\n');
    
    if NA > 1e3
        parallel = 'on'; % PARFOR in TransProbA (Set 'on' if |A| large.)
    else
        parallel = 'off';
    end
    
    P_a = TransProbA( model, parallel );   
    
    fprintf('\t\t\t\t ... DONE!\n');

    %% State space D:
    fprintf('5. Delaunay partition 3D state space as D. ...\n');
    %[D, K] = space_distribution(model);
    
    plot_option = 'none'; 
    warning off all
    [dt, DtriVerts] = simplex2dset_tripart(K, plot_option);
    
    fprintf('\t\t\t\t ... DONE!\n');
    
    %% Intersections of D with all P(a)(Q_k):
    %
    %  * Q_k are partition elements k=1,...,K  of D
    %  * Do this for all action profiles a
    %
    % Usage: 
    %
    % xPolyVerts{a}{k}{i} : action profile a => P(a)
    %                          Current state partition k <=> Q(k)
    %                          i => j(i) <=> Q(j) intersects Q(k)
    % Note: 
    %       * Polygon xPolyVerts{a}{k}{i} in triangle Q(j), 
    %         with number j := xTriIndex{a}{k}(i); and j = 1,...,K.
    %
    %       * Index i = 1,...,N(a,k); where ...
    %         N(a,k) := numel(xPolyVerts{a}{k}) <= K.
    %
    %       * IntersectPolyLcons{a}{k}{i} is linear (weak) inequality
    %         constraints representing polygon xPolyVerts{a}{k}{i}
    
    
    fprintf('6. Index all intersections under P(a)(Q_k) with D. ...\n');
    
    show_waitbar = 'off';
    
    
    tic
        [xQP, xTriIndex, xPolyVerts, xPolyLcons, xPolyRand, xPolyQjRand, ...
            xPolyIndexQjRand] = ...
                    simplex2dset_intersectpmap(NA, P_a, dt, DtriVerts, ...
                                              Nsim, show_waitbar, HitnRun);
    toc

    fprintf('\t\t\t\t ... DONE!\n');
    
    %%
    if strcmp(plot_partsim, 'on')
        % Plot induced mappings: P(a)[Q(j)]:
        plot_options.sort_by_K = 'off';
        plot_options.sort_by_A = 'on';
        simplex2dset_intersectpmap_plot(xPolyVerts, xPolyRand,  ...
                                            xTriIndex, plot_options);
    end
    
    if strcmp(plot_partsiminverse, 'on')
        % Plot Q(j) partitions and random vectors going into each Q(k) that
        % intersects with each P(a)[Q(j)]:
        plot_options.sort_by_K = 'off';
        plot_options.sort_by_A = 'on';
        simplex2dset_intersectpmap_inverse_plot(dt, xPolyQjRand, ...
                                 xTriIndex,xPolyIndexQjRand,plot_options);
    end
    
    %break
    
    %% Main loop: Outer approximation
    fprintf('7. Outer Approximation. ...\n');
    
        % Initial guess: hyperplane levels
        C_old = max(max(H*repmat(V_MAX,np,1)),max(H*repmat(V_MIN,np,1)))...
                                                                *ones(L,K);
        % DIRTY HACK!
        C_old = C_new;
        % Pre-caculate feasible (lambda,b) pairs for every intersection in:
        % { P(a)[Q(j)] \cap Q(k) | a = 1:NA, j,k = 1:K }
        F = feasible_states(model,xPolyQjRand,xPolyIndexQjRand);
    
    %% Set iteration counter for main loop:
    iter_last = 0;				   
	iter = 1;

    %% Initialize epsi dummy:
    epsi = 1;
    
    eps_vec = zeros(1, model.MAXITER1);   % Record distances

    
elseif strcmp(LOAD_OLD,'yes') == 1

    fprintf('MAIN.M: Starting with old saved files!\n');
    % Initial guess of hyperplane levels (payoff correspondence W_outer)
    %load sse_005					% W_old
	iter_last_new = iter_last;
    load(strcat('_mat/sse2014b_iter',int2str(iter_last),'.mat'));
    iter = iter_last_new + 1;
	
end

    %% 
    % Fixed settings for sub GLPK
    N_linear_constraints = L+(NA-1)+1;
    ctype = repmat('U', N_linear_constraints ,1);
    vartype = repmat('C',N_Z,1);
%%
tic


while epsi >= model.TOL1 || iter <= ( iter_last + model.MAXITER1 )

    %% Punishment values $\{ \pi_{k} | k = 1,...,K \}$ given C_old:
    pival = punish_outer( model, xTriIndex, xPolyLcons, ...
                          V_MIN, V_MAX, IndexProfileSet,H, C_old, P_a );

   

    %% Admissibility: W_+ = B(W)
%     C_new = admit_outer( model, xTriIndex, xPolyLcons, V_MIN, V_MAX, ...
%                                 IndexProfileSet, H, C_old, P_a, pival,...
%                                 epsi, iter );
    
                        
 %   tic                    
    C_new = admit_outer_lpset(model, xTriIndex, F, H, C_old, P_a, pival,...
                                              lb, ub, vartype, ...
                                              epsi, iter);
%    toc
                                          
    %% Haussdorff distance:                        
    epsi = mdistance(C_old,C_new);  % d(W(iter),W(iter+1))  
    eps_vec(iter) = epsi;
    
    fprintf('%i\t d(Wo_n,Wo_n+1) = %0.5g\n', iter, epsi);

    %% Intermediate saves for future use / black outs:
    %if epsi < 1
        filename = strcat('_mat/sse2014b_iter',num2str(iter),'.mat');
        save(filename)
    %end

    %% Update guess of C_old:
    C_old = C_new;

    iter = iter + 1;
end

% Calculate time taken
TotalElapsedTime = toc;
AverateTimePerOperator = TotalElapsedTime/(iter - 1 - iter_last);
save _mat/sse2014b_final


% MATLAB profiler option:
if USE_PROFILER == 1
    profile viewer
    p = profile('info');
    profsave(p,'_ProfileResults')
end
