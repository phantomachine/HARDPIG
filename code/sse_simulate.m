%% SSE_SIMULATE.M
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

 % Add path to directory:
addpath(strcat(pwd,'/glpkmex/'));
addpath(strcat(pwd,'/_func/'));
addpath(strcat(pwd,'/_func/sphereFit'));
addpath(strcat(pwd,'/eq_sphere_partitions/eq_partitions/'));
addpath(strcat(pwd,'/compecon2011/CEtools/'));
addpath(strcat(pwd,'/_mat/'));

load_stratpoly = 1;     % Load save PolyIntersections

T = 30;
Nagrid = 5;
Nbgrid = 6;

smoother = 0.1;        % Smoother if stuck

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
lambda_now = [0.2081    0.0010    0.7909];

%% INNER APPROX & CENTROIDS OF POLYTOPE SLICES FOR EACH Q(k)
plothull = 'off';
stepcorr = stepcorrchull(H, C_new, K, plothull);

VertSetChull = stepcorr.VertSetChull;
G_normal = stepcorr.Gsearch;
c_normal = stepcorr.c_normal;
centroids = stepcorr.Wcentroid_hull;

%% INNER RAY: FINER SEARCH SUBGRADIENTS
% Fix search subgradients; Use Leopardi (2006):
M_theta = 400;                              % Inner ray subgradients
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
    fprintf('\n Initial  w_now is inside W(k), k = %0.5g\n',k_now);
else
    error('STRAT_INNERAY:IsInside','Initial  w_now NOT inside W(k)')
end

%% SIMULATE SSEs
for t = 1:T
    %% Feasible set of policies b, B(lambda)
    Feas_temp = lambda_now * model.ProfileSetB';
    %B_temp = ProfileSetB( (Feas_temp >= -tolB & Feas_temp <= tolB), :)
    B_temp = model.ProfileSetB( (Feas_temp <= tolB), :);
    NB = size(B_temp,1);
    
    %% (1) Check if w_now is an extreme point supported by pure strategies
    % Evaluate R() function in ray direction wnow_ray:
    wko = centroids(k_now,:)'; % Column vector
    wnow_ray = (w_now - wko)/norm(w_now - wko,2);   
    [l_star, aopt_idx, bopt_idx,~] = pureray(model, wnow_ray, ...
                           lambda_now, k_now, B_temp, P_a, ...
                           wko, dt, DtriVerts, G_normal, c_normal, pival );
                       
    fprintf('\n t = %i \tCheck W_NOW extreme point ... l* = %6.4d',t,l_star);
    
    if l_star == 1.0
        fprintf('\n t = %i \t OK! w_now is an extreme point.', t);
        %% CASE 1: w_now is an extreme point itself
        % ----------------------------------------------------------------
        a_star = model.ProfileSetA(aopt_idx,:);
        b_star = B_temp(bopt_idx,:);
        
        theta_star = wnow_ray;
        %         % Current payoff from a: effort disutility
        %             Util_a = -phia(model, a_star);
        %         % Current payoff from b:
        %             cons = b_star;
        %             cons(N+1:end) = cons(N+1:end) + model.WAGE;
        %             Util_b = ucons(model, cons);
        %         % Total current payoff: (column vector!)
        %             Util_now = model.DELSTAR*( Util_a + Util_b )';
        %         % Continuation promised values:
        %         w_next = (wko + l_star*wnow_ray - Util_now)/model.DELTA;
    else
        fprintf('\n\tNow finding extreme points to randomize over ...');
        %% CASE 2: Support w_now by randomizing over all extreme points
        % ----------------------------------------------------------------
        % (2) We need to randomize over extreme points containing w_now:
        lpos_temp = nan(M_theta,1);
        tpos_temp = nan(M_theta,1);
        apos_temp = nan(M_theta, N_Z);
        bpos_temp = nan(M_theta, N_Z);
        reverseStr = '';
        for theta_idx = 1:M_theta
                        
            % Current ray direction
            theta_ray = Theta(theta_idx,:)';
            % Check if exist, then store extreme points:
            [l_star,aopt_idx,bopt_idx,L_m]= pureray(model, theta_ray, ...
                               lambda_now, k_now, B_temp, P_a, wko, dt, ...
                                     DtriVerts, G_normal, c_normal, pival);            
            
            if (l_star > 0.0) && (l_star ~= inf)
                % Store (a,b,theta,l) in temp storage:
                tpos_temp(theta_idx) = theta_idx;
                lpos_temp(theta_idx) = l_star;
                apos_temp(theta_idx,:) = model.ProfileSetA(aopt_idx,:);
                bpos_temp(theta_idx,:) = B_temp(bopt_idx,:);           
            end
            
            %fprintf('\n\tEvaluating PURERAY at theta_idx = %i of %i', ...
            %                                          theta_idx, M_theta);
            % Display the progress
            percentDone = 100 * theta_idx / M_theta;
            msg = sprintf('\t %3.1f Percent Done ...', percentDone); 
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));

        end
        
        % Extreme Points and supporting pure strategies:
        tpos = tpos_temp(~isnan(tpos_temp)); % Indices to Theta elements
        M_ext = length(tpos); % No. of extreme points (vertices)
        
        if isempty(tpos)
            error('tpos is empty! There are no solutions')
        else
            fprintf('\n\tFinding probs.over M* = %i vertices',M_ext);
        end
        
        lpos = lpos_temp(tpos); % l* supporting vertices
        
        apos = apos_temp(tpos, :); % a* (pure actions)
        bpos = bpos_temp(tpos, :); % b* (pure actions)
        
        theta_ext = Theta(tpos, :); % rays/directions of vertices
        Z_ext = (repmat(wko', M_ext, 1) + repmat(lpos,1,N_Z).* theta_ext)'; 
                                                                % vertices
        
        % (3) Now given stored extreme points, find probability weights:
            Objectif = ones(1,M_ext);
            Constraints_lhs = [ Z_ext;             % 'U' -- inequalities
                                ones(1,M_ext) ];   % 'S' == equality
            Constraints_rhs = [ w_now;
                                1             ];
            lb = zeros(M_ext,1);
            ub = ones(M_ext,1);
            ctype = [ repmat('U', N_Z,1); 'S' ];
            vartype = repmat('C', M_ext,1);
            sense = -1;
            param.msglev=1; % Output GLPK messages
            param.save=1;   % Set save options
            [probs, ~, ~, ~] = glpk (Objectif, Constraints_lhs, ...
                                Constraints_rhs, lb, ub, ctype, vartype,...
                                sense, param);
            selection = (probs ~= 0); % Select vertices with +ve weights
            probs = probs(selection);
            apos = apos(selection, :);
            bpos = bpos(selection, :);
            lpos = lpos(selection, :);
            theta_ext = theta_ext(selection,:);
            
            csumprob = cumsum(probs);
            
        % (4) Generate a uniform r.v. "sunspot" realization
            draw = rand(1);
            % Then see where DRAW falls in the support of CSUMPROB:
            for i = 1:length(csumprob)-1
                if i == 1 && draw <= csumprob(1)
                    x = i;
                elseif draw > csumprob(i) && draw <= csumprob(i+1)
                    x = i+1;
                end
            end
            
        % (5) Given sunspot x, pick pure actions and continuation payoffs:
            a_star = apos(x,:);
            b_star = bpos(x,:);
            l_star = lpos(x);
            theta_star = theta_ext(x,:)';
            
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
    fprintf('\n\t\t\t\tUpdating actions and continuation values ...');
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
