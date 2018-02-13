%% STRAT_SUNSPOT.M
%
% Construct sample pure strategies supporting initial (lambda0, v0) given
% inner approximating correspondence W^i := [ G(l,k), c(l,k) ]. Uses Bilinear
% Programming approach.
%
% (c) 2012, 2013, Timothy Kam. Email: tcy.kam__at__gmail.com 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Use subject to GNU LGPL licensing terms. Cite this header and 
% author completely in subsequent re-use and modifications.
% =========================================================================
% $Revision: 5.0.0 $  $Date: 2013/09/11 12:45:20 $

close all
clear all
clc

UNIX = 'on';

 % Add path to directory:
    if strcmp(UNIX,'on') == 1   
        addpath(strcat(pwd,'/glpkmex/'));
        addpath(strcat(pwd,'/_func/'));
        addpath(strcat(pwd,'/eq_sphere_partitions/'));
        addpath(strcat(pwd,'/_mat/'));
    else
        addpath(strcat(pwd,'\glpkmex\'));
        addpath(strcat(pwd,'\_func\'));
        addpath(strcat(pwd,'\eq_sphere_partitions\'));
        addpath(strcat(pwd,'\_mat\'));
    end
    
%% SETTINGS
load_old = 0;

    load_stratpoly = 1;     % Load save PolyIntersections
    
    T = 10;
    Nagrid = 5;
    Nbgrid = 6;
    
    smoother = 0.9;        % Smoother if stuck

    % Small tolerance for GBC deficit/surplus
    tolB = 0.15;

    % GLPK options
    param.msglev=0; % Output GLPK messages on workspace      
    param.save=0;   % Set save options


%% SAVED RESULTS    
 
    % Load saved results from SSE
    load('_mat/sse2014_iter162.mat')
    
    % load Fixed policy result
    load('_mat/ramsey_lr01.mat')
    
    %% Initial state of interest
    lambda_now = lambda_ramsey_ss_max;
    
    %% CENTROIDS OF POLYTOPE SLICES FOR EACH Q(k), k = 1,...,K
    plothull = 'off';
    stepcorr = stepcorrchull(H, C_new, K, plothull);

    VertSetChull = stepcorr.VertSetChull;
    Gsearch = stepcorr.Gsearch;
    c_normal = stepcorr.c_normal;
    
%% NEW SHIT ....

%% RICHER ACTION PROFILE SETS FOR SIMULATING STRATEGIES
if strcmp(load_old, 'on') == 0
    
    [ ProfileSetA, ProfileSetB, NA, NB ] = strat_abspace(model,...
                                               Nagrid, Nbgrid, N_Z,M,N);

    model.ProfileSetB = ProfileSetB;
    model.ProfileSetA = ProfileSetA;
 
%% REDO TRANSITION PROBABILITIES WITH RICHER SET A
    parallel = 'on';
    P_a = TransProbA(model,parallel);
    
    DELSTAR = (1-DELTA);

%% REDO STATE SPACE INTERSECTIONS
   
    warning off all
    if load_stratpoly == 1
        load('_mat/stratpoly.mat')
    else

        fprintf('Index all intersections under P(a)(Q_k) with D. ...\n');

            show_waitbar = 'off';


            tic
                [xQP, xTriIndex, xPolyVerts, xPolyLcons, xPolyRand, ...
                    xPolyQjRand, xPolyIndexQjRand] = ...
                     simplex2dset_intersectpmap(NA, P_a, dt, DtriVerts, ...
                                       Nsim, show_waitbar, HitnRun);
            toc

        fprintf('\t\t\t\t ... DONE!\n');
    end


%% CONSTRUCT SAMPLE SSE PATH --- INDEXES TO (a,b) in A x B

    % Storage
    alpha = cell(T,1);  % Strategies a
    beta = cell(T,1);   % Strategies b
    wvec = beta;        % Promised Values
    lams = beta;        % Distributions

    val_a = zeros(1,NA);    % Temp storage for loop a_idx
    wnext_a = cell(NA,1);
    lnext_a = wnext_a;
                              
    % Initial conditions
    %lambda_now = [1/3,1/3,1/3];
    
    k_now = InPartElement(lambda_now, dt, DtriVerts); % initial state index
    
    % Convex hull of extreme points of set W(:,k_now)
    X = VertSetChull{k_now};    % Set of extreme points of W(k_now)
    
    
    %% Maximal Government payoff given (lambda_now, X)
    [ vG_max, index_max ]  = max(lambda_now * X');
    w_now = X(index_max,:)';
    
        % Check w_now is contained in W(k_now)
        tri = delaunayn(X);         % 3D-Delaunay triangulations
        tn = tsearchn(X, tri, w_now'); 
        IsInside = ~isnan(tn);
        
        if IsInside
            sprintf('\n Initial  w_now is inside W(k), k = %0.5g\n',k_now);
        else
            error('STRAT_INNERAY:IsInside',...
                   'Initial  w_now NOT inside W(k)')
        end
        

%% GENERATE SUNSPOTS
seed = 1;

extra = ceil(T/2);

[ x, seed ] = uniform_on_simplex01_map ( L, T + extra, seed );
x = x';

%% MAIN LOOP: COMPUTING STRATEGIES SUPPORTING (lambda_now, w_now)   


% % %% HACK FOR NOW!
% % lambda_now = [0.2081    0.0010    0.7909];
% % w_now = [ -0.6096;     1.7172;     0.8388 ];

% Sunspot realization at date t = 1   

fprintf('Simulating Strategies by Sunspots ... please wait\n');

t = 1;
xspot = x(t,:);
    
while t <= T
    
    % Current utilitarian social payoff
    v_now = lambda_now * w_now;
   
    %% Work out (approximately) feasible set of policies b, B(lambda)
    
    Feas_temp = lambda_now * ProfileSetB';
    B_temp = ProfileSetB( (Feas_temp <= tolB), :);
    NB = size(B_temp,1);
    
    % Variable storage:
    val_b = zeros(NB,1);
    aoptindex_b = cell(NB,1);
    wnext_b = cell(NB,1);
    lnext_b = cell(NB,1);
        
    % Inner loops:
    for b_idx = 1 : NB             % Loop of feasible policies b in B_temp
           
        parfor a_idx = 1 : NA       % loop over agents' actions a in A
            
            xspot_local = xspot;
            
            lb_local = lb;
            ub_local = ub;
            param_local = param;
            vartype_local = vartype;
            NZ_local = N_Z;
            
            ProfileSetA_local = ProfileSetA;
            P_a_local = P_a;
            B_temp_local = B_temp;
            
            H_local = H;                % Search subgradients
            %c_new_local = C_new;
            
            G_normal_local = Gsearch;   % Approx subgradients (normals)
            c_normal_local = c_normal;
            
            delta_local = DELTA;
            delstar_local = DELSTAR;
            
            v_now_local = v_now;
            Vmin = V_MIN;
            Vmax = V_MAX;
            
            % Current payoff from a: effort disutility
            Util_a = -phia(model,ProfileSetA(a_idx,:)); 
            Util_ia = repmat(Util_a, NA-1, 1);
            
            % All other a# \neq a, and induced P(a#) matrices:
            Index_A = 1:NA;                     % indexes to { a }
            ia_Not = setdiff(Index_A, a_idx); 
                                                % All a# \neq a
            P_ia_Not = P_a_local(:,:,ia_Not);   % All P(a#) \neq P(a), 
                                                % N x N x card(A)-1
            P_ia_not = reshape(P_ia_Not,...
                                N_Z*(NA-1),N_Z);
                                            % Tile this too!
                                            
            % Current payoffs from all deviations from a:        
            Util_ia_Not = -phia(model, ProfileSetA_local(ia_Not,:) );
            
            % Transitions
            P_now = P_a_local(:,:,a_idx);
            P_ia = repmat(P_now,NA-1,1);
            
            % Agents' incentive constraints at current {k,l,a}
            DevP = delta_local*(P_ia_not - P_ia);
            DevA = delstar_local*(Util_ia - Util_ia_Not);
            [nr, nc ] = size(DevA);
            DevA = reshape(DevA', nr*nc, 1);
            
            % Current payoff from b:
            b_current_feasible = B_temp_local(b_idx,:);
            cons = b_current_feasible;
            cons(N+1:end) = cons(N+1:end) + model.WAGE;
            Util_b = ucons(model, cons);  

            % Total current payoff:        
            Util_now = delstar_local*( Util_a + Util_b )'; 
                                            % column vector
                           
          
                           
           % Check which Q(k) contains lambda_now:
           %k_now = InPartElement(lambda_now, dt, DtriVerts);
                           
           % Get continuation game state:
            
            lambda_next = lambda_now * P_now;
            
            % Check which Q(k) contains lambda_next:        
            k_next = InPartElement(lambda_next, dt, DtriVerts);
            
            % Get appropriate partitition element's correspondence:
            %c_next = C_new_local(:,k_next);
            G_next = G_normal_local{k_next};
            c_next = c_normal_local{k_next};
            %clevel = C_new_local(:,k_next);
            
            % Solve LP problem for maximizers w_next that enforces current
            % action-policy pair (a,b):
            %
            % NOTE: (1)-(2): weak inequality constraints
            %       (3) is strict equality constraint
            
            constraint = [  G_next;                                %... (1)
                            DevP;                                  %... (2)
                            lambda_now * delta_local*P_now   ];    %... (3)
                                        
            bconstraint1 = [ c_next;                              %...  (1)
                            DevA                      ];          %...  (2)
            bconstraint2 =  v_now - lambda_now * Util_now ;       %...  (3)
                                          
                bconstraint = [ bconstraint1; bconstraint2 ] ;
                
            % Integral over sunspot xspot := weighted avg. over vertices
            %
            %           (1 x L)     (L x N_Z)               (N_Z x N_Z)
            objective = xspot_local *H_local   *delta_local    *P_now;
            
            %lambda_next; %ones(1,N_Z);
            
                    % Constratint types
                    ctype_local1 = repmat('U', size(bconstraint1,1),1);
                                                        % Inequality
                                                        % constraints
                    ctype_local2 = repmat('S', size(bconstraint2,1),1);
                                                        % Equality
                                                        % constraints
                    
                    ctype_local = [ ctype_local1; ctype_local2 ];
                    
                    %vartype = repmat('C',N_Z,1);

                    %% GLPK                       
                    %param.msglev=1; % Output GLPK messages     
                    %param.save=1;   % Set save options
                    sense = -1;     % Max (-1) or Min (+1) problem

                    [w_temp,fval,status,~] = glpk( objective,...
                                                constraint,...
                                                bconstraint,...
                                                lb_local,ub_local,...
                                                ctype_local,...
                                                vartype_local,...
                                                sense,param_local);
                    
                    wnext_a{a_idx} = w_temp;
                    lnext_a{a_idx} = lambda_next;
                    
                    if status == 2     
                        %status_vec(n) = status; 
                        %warning('STRAT:status',...
                        %            '[2] Feasible but not optimal');
                        val_a(a_idx) = xspot_local*H_local*Util_now + fval;
                    elseif status == 5
                        %status_vec(n) = status;
                        val_a(a_idx) = xspot_local*H_local*Util_now + fval;
                    else
                        val_a(a_idx) = -inf;
                        %warning('STRAT:status',...
                        %                      'Solution not found');
                    end          
                           
        end % EndFor a
                
        % Given b, maximize over the a's:
        
        [ vmax_temp, ind_temp  ] = max(val_a);
        
        if isinf(vmax_temp)
            warning('STRAT:vmax_temp','No feasible solution found')
        end
            val_b(b_idx) = vmax_temp; 
            aoptindex_b{b_idx} = ind_temp;

            wnext_b{b_idx} = wnext_a{ind_temp};
            lnext_b{b_idx} = lnext_a{ind_temp};
        
        %fprintf('t = %i \t b_idx = %i \n',t, b_idx);

        
    end % EndFor b
    
    % Maximize over the feasible b's:
    
    [vtemp, boptindex ] = max(val_b);
    
    fprintf('t = %i \t vG = %0.5g \t bopt_idx = %i \t aopt_idx = %i \n',...
                             t, vtemp, boptindex, aoptindex_b{boptindex});
    
    if isinf(vtemp)
        warning('STRAT:boptindex','No feasible solution found')
        
        % Then redo stage t with new sunspot randomization
        
        if t+1 <= T + extra
            xspot = x(t+1,:);
        else
            warning('STRAT:xspot','You have run out of sunspots!')
            warning('STRAT:xspot','... so I am recycling the first one!')
            xspot = x(1,:);
        end
        
        % Re-Initialize conditions by small perturbation
        k_now = InPartElement(lambda_now, dt, DtriVerts); 
                                                    % initial state index
    
        % Convex hull of extreme points of set W(:,k_now)
        hullvert = VertSetChull{k_now};    % Set of extreme points of W(k_now)
    
    
        %% Maximal Government payoff given (lambda_now, hullvert)
        [ ~, index_max ]  = max(lambda_now * hullvert');
        z = hullvert(index_max,:)';
        
        %% Take a mixture of w_now and extreme point z
        w_now = smoother*z + (1-smoother)*w_now;
        
        % And now t is still t, so we repeat!
        
    else
        
        % Store successful solutions:
        
        w_next = wnext_b{boptindex};
        %l_opt = lnext_b{boptindex};

        % Store profiles (a*, b*) at stage/date t:    
        alpha{t} = ProfileSetA(aoptindex_b{boptindex},:); 
                                                 % Index to row vector in A
        beta{t} = B_temp(boptindex,:);        
                                            % Index to row vector in B_temp

        % Update state and continuation payoffs
        lambda_next = lambda_now * P_a(:,:,aoptindex_b{boptindex});

        lambda_now = lambda_next;
        w_now = w_next;

        wvec{t} = w_now;
        lams{t} = lambda_now;
        
        % Update time counter
        t = t + 1;
        
        % Sunspot realization at date t    
        xspot = x(t,:);
    end
        
    
    
end % EndFor t


%% Pack it all in:
strat.a = alpha;
strat.b = beta;
strat.w = wvec;
strat.l = lams;
                           
save(strcat('_mat/stratsunspot',int2str(T),'.mat'), 'strat')

else

load(strcat('_mat/stratsunspot',int2str(T),'.mat'), 'strat')
 
    sort_ascend_N = (N:-1:1);
    
    beta = cell2mat(strat.b);
        beta(:,1:N) = beta(:, sort_ascend_N );
    alpha = cell2mat(strat.a);
        alpha(:,1:N) = alpha(:, sort_ascend_N );
    distros = cell2mat(strat.l);
        distros(:,1:N) = distros(:, sort_ascend_N );
    values = cell2mat(strat.w);
        values = reshape(values, T, N_Z);
        values(:,1:N) = values(:, sort_ascend_N );
    
    % Plot data:
    
    x = [-N:-1,1:M];
    y = 1:T;
    
    z = cell(4,1);
    z{1} = beta;
    z{2} = alpha;
    z{3} = distros;
    z{4} = values;
    
    %xy_zero = gridmake([-N; M],[1; T], 0)';
    x_zero = [ -N, -N, M, M ];
    y_zero = [  1,  T, T, 1 ];
    z_zero = zeros(1,4);
 
    znames = { 'Transfer', 'Effort', 'Proportion', 'ADP' };
    
    rgb = [0.6 0.8 1.0];
    
    for i = 1:numel(z)
        figure
            fill3(x_zero,y_zero,z_zero,rgb,'FaceAlpha', 0.2,...
                                                    'edgecolor','none')

            hold on
                stem3(x,y,z{i}, ':.b')
                set(gca,'XTick',x)
            hold off
            %set(gca,'XTickLabel',['-2';'-1';'+1'])
            box off
            grid off
            xlabel('Z')
            ylabel('Time')
            zlabel(znames{i})
            
            fname = strcat('_figures/','strat-',znames{i});
            hgsave(fname)
            print('-depsc',fname)
    end

end

%% TABULATE RESULTS

% Select time periods:
t = [ 3, 10, 25, 28, 30 ];

   matrix = [t', beta(t,:), alpha(t,:), distros(t,:)];
   
   rowLabels = {'t'};
   columnLabels = {'$t$', '$b(-2)$', '$b(-1)$', '$b(+1)$',...
                            '$a(-2)$', '$a(-1)$', '$a(+1)$',...
                        '$\lambda(-2)$', '$\lambda(-1)$', '$\lambda(+1)$'};
   matrix2latex(matrix, '_figures/strat-table.tex', ...
                        'columnLabels', columnLabels, ...
                        'alignment', 'c', ...
                        'format', '%-0.2g', 'size','small');
%
    