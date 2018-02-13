%% SSE_SIMULATE_PARALLEL.M
%
% Construct sample pure strategies supporting initial (lambda0, v0) given
% inner approximating correspondence W^i := [ G(l,k), c(l,k) ].
% Uses GLPK-MEX interface by Nicolo' Giorgetti.
%
% (c) 2012, 2013, Timothy Kam. Email: tcy.kam__at__gmail.com
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Use subject to GNU LGPL licensing terms. Cite this header and
% author completely in subsequent re-use and modifications.
% =========================================================================
% $Revision: 5.0.0 $  $Date: 2015/12/24 23:45:20 $

close all
clearvars
clc
warning off all


%% SETTINGS
load_old = 'off';
more_action = 'off';
plot_figures = 'on';
enrich_thetaset = 'off';

% Add path to directory:
addpath(strcat(pwd,'/glpkmex/'));
addpath(strcat(pwd,'/_func/'));
addpath(strcat(pwd,'/_func/sphereFit'));
addpath(strcat(pwd,'/eq_sphere_partitions/eq_partitions/'));
addpath(strcat(pwd,'/compecon2011/CEtools/'));
addpath(strcat(pwd,'/_mat/'));

% Simulation length (Time)
T = 20;

% Option to redefine gridpoints in action-policy (a,b) spaces
Nagrid = 2;
Nbgrid = 2;

% Small tolerance tolB% for GBC deficit/surplus
tolB = 0.15;

% GLPK options
param.msglev=0; % Output GLPK messages on workspace
param.save=0;   % Set save options

%% SAVED RESULTS

% load Fixed policy result
load('_mat/ramsey_lr01.mat')

% Load saved results from SSE
load('_mat/sse2014_iter162.mat')



%% Initial state of interest
%lambda_now = lambda_ramsey_ss_max;
%lambda_now = [0.2081    0.0010    0.7909];
%lambda_now = [0.099, 0.001, .9];
%lambda_now = [.108, .3, .592];
lambda_now = [.1, .3, .6];

if strcmp(load_old, 'on')
    load(strcat('_mat/sse_simulate',int2str(T),'.mat'), 'strat')
else
    if strcmp(more_action, 'on')
        % Richer action-policy spaces
        [ ProfileSetA, ProfileSetB, NA, NB ] = strat_abspace(model,...
                                                 Nagrid, Nbgrid, N_Z,M,N);

        model.ProfileSetB = ProfileSetB;
        model.ProfileSetA = ProfileSetA;

        % Redo transitions
        parallel = 'on';
        P_a = TransProbA(model,parallel);
    end

    %% INNER APPROX & CENTROIDS OF POLYTOPE SLICES FOR EACH Q(k)
    plothull = 'off';
    stepcorr = stepcorrchull(H, C_new, K, plothull);

    VertSetChull = stepcorr.VertSetChull;
    G_normal = stepcorr.Gsearch;
    c_normal = stepcorr.c_normal;
    centroids = stepcorr.Wcentroid_hull;

    %% INNER RAY: FINER SEARCH SUBGRADIENTS
    % Fix search subgradients; Use Leopardi (2006):
    M_theta = 1500;                              % Inner ray subgradients
    Theta = subgradient(np,M_theta);            % H is finite subset of the
                                                % np-dimensional sphere:
                                                % S^(np-1). So, for each i,
                                                % sum( H(i,:).^2, 2 ) = 1.

    %% CONSTRUCT SAMPLE SSE PATH --- INDEXES TO (a,b) in A x B
    % Storage
    alpha = cell(T,1);  % Strategies a
    beta = cell(T,1);   % Strategies b
    wvec = cell(T+1,1);   % Promised Values
    lams = cell(T+1,1);   % Distributions

    % map lambda_now into k \in Q_k
    k_now = InPartElement(lambda_now, dt, DtriVerts); % initial state index

    % Convex hull of extreme points of set W(:,k_now)
    X = VertSetChull{k_now};    % Set of extreme points of W(k_now)

    % Maximal Government payoff given (lambda_now, X)
    [ vG_max, index_max ]  = max(lambda_now * X');
    w_now = X(index_max,:)';
    wvec{1} = w_now;
    lams{1} = lambda_now;

    % Check w_now is contained in W(k_now)
    tri = delaunayn(X);         % 3D-Delaunay triangulations
    tn = tsearchn(X, tri, w_now');
    IsInside = ~isnan(tn);

    if IsInside
        fprintf('\n Initial  w_now is inside W(k), k = %i\n',k_now);
    else
        error('STRAT_INNERAY:IsInside','Initial  w_now NOT inside W(k)')
    end

    %parpool

    %% SIMULATE SSEs
    for t = 1:T
        %% Feasible set of policies b, B(lambda)
        Feas_temp = lambda_now * model.ProfileSetB';
        B_temp = model.ProfileSetB( (Feas_temp >= -tolB | Feas_temp <= tolB), :);
        %B_temp = model.ProfileSetB( (Feas_temp <= tolB), :);
        NB = size(B_temp,1);

        %% (1) Check if w_now an extreme point supported by pure strategies
        % Evaluate R() function in ray direction wnow_ray:
        wko = centroids(k_now,:)'; % Column vector
        wnow_ray = (w_now - wko)/norm(w_now - wko,2);
        [l_star, aopt_idx, bopt_idx,~] = pureray(model, wnow_ray, ...
                               lambda_now, k_now, B_temp, P_a, ...
                          wko, dt, DtriVerts, G_normal, c_normal, pival );

        fprintf('\n t = %i | k = %i | l* = %6.4g',t, k_now, l_star);

        if l_star == 1.0
            fprintf('\n t = %i \t OK! w_now is an extreme point.', t);
            %% CASE 1: w_now is an extreme point itself
            % -------------------------------------------------------------
            a_star = model.ProfileSetA(aopt_idx,:);
            b_star = B_temp(bopt_idx,:);

            theta_star = wnow_ray;
        else
            fprintf('\n\t\tFinding extreme points to randomize over ...');
            %% CASE 2: Support w_now by randomizing over all extreme points
            % -------------------------------------------------------------
            % (2) We need to randomize over extreme points containing w_now
            
            %reverseStr = '';
            
            % Enrich Theta with G_normal{k}:
            if strcmp(enrich_thetaset, 'on')
                if length(H) >= length(G_normal{k_now})
                    G_temp = [ H; setdiff(H, G_normal{k_now}, 'rows')];
                else
                    G_temp = [ G_normal{k_now}; ...
                               setdiff(G_normal{k_now},H,'rows')];
                end

                if length(Theta) >= length(G_temp)
                    H_temp = [ Theta; setdiff(Theta, G_temp, 'rows')];
                else
                    H_temp = [ G_temp; ...
                               setdiff(G_temp, Theta,'rows')];
                end
                %H_temp = setdiff(Theta, G_temp, 'rows');
                Theta_now = [Theta; H_temp];
            else
                Theta_now = Theta;
            end
            
            M_theta_now = size(Theta_now, 1);
            
            lpos_temp = nan(M_theta_now,1);
            tpos_temp = nan(M_theta_now,1);
            apos_temp = nan(M_theta_now, N_Z);
            bpos_temp = nan(M_theta_now, N_Z);
            
            parfor theta_idx = 1:M_theta_now
                %Theta_local = Theta_now;
                % Current ray direction
                theta_ray = Theta_now(theta_idx,:)';
                % Check if exist, then store extreme points:
                [l_star,aopt_idx,bopt_idx,L_m]= pureray(model,theta_ray,...
                               lambda_now, k_now, B_temp, P_a, wko, dt, ...
                                     DtriVerts, G_normal, c_normal, pival);

                if (l_star > 0.0) && (l_star ~= inf)
                    % Store (a,b,theta,l) in temp storage:
                    tpos_temp(theta_idx) = theta_idx;
                    lpos_temp(theta_idx) = l_star;
                    apos_temp(theta_idx,:) = model.ProfileSetA(aopt_idx,:);
                    bpos_temp(theta_idx,:) = B_temp(bopt_idx,:);
                end
            end

            % Extreme Points and supporting pure strategies:
            tpos = tpos_temp(~isnan(tpos_temp)); %Indices to Theta elements
            M_ext = length(tpos); % No. of extreme points (vertices)

            if isempty(tpos)
                error('tpos is empty! There are no solutions')
            else
                fprintf('\n\tFinding probs.over M* = %i vertices',M_ext);
            end

            lpos = lpos_temp(tpos); % l* supporting vertices

            apos = apos_temp(tpos, :); % a* (pure actions)
            bpos = bpos_temp(tpos, :); % b* (pure actions)

            theta_ext = Theta_now(tpos, :); % rays/directions of vertices
            Z_ext = (repmat(wko',M_ext,1)+repmat(lpos,1,N_Z).*theta_ext)';
                                                                % vertices

            % (3) Now given stored extreme points, find probability weights
                Objectif = zeros(1,M_ext);
                Constraints_lhs = [ Z_ext;             % 'U' --inequalities
                                    ones(1,M_ext) ];   % 'S' == equality
                Constraints_rhs = [ w_now;
                                    1             ];
                lb = zeros(M_ext,1);
                ub = ones(M_ext,1);
                ctype = [ repmat('S', N_Z,1); 'S' ];
                vartype = repmat('C', M_ext,1);
                sense = 1;
                param.msglev=0; % Output GLPK messages
                param.save=0;   % Set save options
                [probs, ~, status, extra] = glpk(Objectif, Constraints_lhs, ...
                                    Constraints_rhs,lb,ub,ctype,vartype,...
                                    sense, param);
                selection = (probs ~= 0); %Pick vertices w/ +ve weight
                if ~isempty(selection)
                    probs = probs(selection);
                    fprintf('\n\t\t\tRandomizing over %i vertices', ...
                                                                numel(probs));
                    apos = apos(selection, :);
                    bpos = bpos(selection, :);
                    lpos = lpos(selection, :);
                    theta_ext = theta_ext(selection,:);
                else % Uniform weights over all extreme points (ad-hoc)
                    probs = ones(M_ext,1)/M_ext;
                    fprintf('\n\t\t\tRandomizing over all M* = %i vertices',...
                                                                        M_ext);
                end
                csumprob = cumsum(probs);

                if numel(csumprob) == 1         % degenerate lottery
                    % (4) No need to randomize
                    fprintf('\n\t\t\tDegenerate lottery on M* = %i vertices',...
                                                                        M_ext);
                    a_star = apos;
                    b_star = bpos;
                    l_star = lpos;
                    theta_star = theta_ext';
                else % nondegenerate lottery
                    % (4) Generate a uniform r.v. "sunspot" realization
                    draw = rand(1);
                    fprintf('\n\t\t\t\tSunspot draw = %4.2g', draw);
                    % Then see where DRAW falls in the support of CSUMPROB:
                    for i = 1:length(csumprob)-1
                        if i == 1 && draw <= csumprob(1)
                            x = i;
                        elseif draw>csumprob(i) && draw<=csumprob(i+1)
                            x = i+1;
                        end
                    end
                    % (5) Given sunspot x, pick pure actions and
                    % continuation payoffs:
                    a_star = apos(x,:);
                    b_star = bpos(x,:);
                    l_star = lpos(x);
                    theta_star = theta_ext(x,:)';
                end
        end
        % Current payoff from a: effort disutility
        Util_a = -phia(model, a_star);
        % Current payoff from b:
        cons = b_star;
        cons(N+1:end) = cons(N+1:end) + model.WAGE;
        Util_b = ucons(model, cons);
        % Total current payoff: (column vector!)
        Util_now = model.DELSTAR*( Util_a + Util_b )';
        w_next = (wko + l_star*theta_star - Util_now)/model.DELTA;

        %% STORE RESULTS: (a,b,w+,lambda+)
        fprintf('\n\tUpdating actions and continuation values ...\n');
        alpha{t} = a_star;
        beta{t} = b_star;
        wvec{t+1} = w_next;
            w_now = w_next; % Update to start next period
        astar_idx = getindex(a_star, model.ProfileSetA);
        lambda_next = lambda_now * P_a(:,:,astar_idx);
        lams{t+1} = lambda_next;
            lambda_now = lambda_next; %Update to start next period
            k_now = InPartElement(lambda_now, dt, DtriVerts); % Q_k
    end

    delete(gcp) % Close parpool

    %% Pack it all in:
    strat.a = alpha;
    strat.b = beta;
    strat.w = wvec;
    strat.l = lams;

    save(strcat('_mat/sse_simulate',int2str(T),'.mat'), 'strat')

end

%% PLOT RESULTS

if strcmp(plot_figures, 'on')
    beta = cell2mat(strat.b);
    alpha = cell2mat(strat.a);
    distros = cell2mat(strat.l);
    distros = distros(1:T,:);
    values = cell2mat(strat.w);
    values = reshape(values, T+1, N_Z);
    values = values(1:T,:);

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

                fname = strcat('_figures/','sse_simulate-',znames{i});
                hgsave(fname)
                print('-depsc',fname)
        end
end
