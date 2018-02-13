%% STRAT_INNERRAY.M
%
% Construct sample pure strategies supporting initial (lambda0, w0) given
% inner approximating correspondence W^i := [ G(l,k), c(l,k) ]. This adapts
% the inner-ray method proposed by Judd, Yeltekin and Conklin (2003, Ecta).
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

%% SETTINGS

T = 4;

Nagrid = 10;
Nbgrid = 10;

% Small tolerance for GBC deficit/surplus
tolB = 0.15;

% GLPK options
    param.msglev=0; % Output GLPK messages on workspace      
    param.save=0;   % Set save options

%% SAVED RESULTS    
    
% load Fixed policy result
load('_mat/ramsey_lr01.mat')

lambda_now = lambda_ramsey_ss_max;

wrss_now = vrmax;


% Load saved results from SSE
load('_mat/sse2014_iter162.mat')

    L = size(lambda_new, 3);
    K = size(lambda_new, 4);

    % lambda_new and w_new are (1 x 3 x L x K ) arrays
    vG_k = sum( lambda_new .* w_new, 2 );
    vG_k = reshape(vG_k, L, K);
    [vG_k_max, max_index] = max(vG_k,[],1);
    [vG_k_min, min_index] = min(vG_k,[],1);

    wmin = zeros(K,N_Z);
    wmax = wmin;
    lmin = wmax;
    lmax = lmin;

    for k = 1 : K
        wmin(k,:) = w_new(:,:,min_index(k),k)'; % w supporting min(vG)
        wmax(k,:) = w_new(:,:,max_index(k),k)'; % w supporting max(vG)
        lmin(k,:) = lambda_new(:,:,min_index(k),k)';
                        % Not needed: lamda supporting min(vG)
        lmax(k,:) = lambda_new(:,:,max_index(k),k)';
                        % Not needed: lamda supporting max(vG)
    end


%% CENTROIDS OF POLYTOPE SLICES FOR EACH Q(k), k = 1,...,K

    % Vertices of W(k), fit approximate centroid:
    VertSet = cell(K,1);
    VertSetChull = VertSet;
    Vert_chull_ind = VertSet;
    Gsearch = VertSet;
    c_normal = VertSet; 
    Wcentroid = zeros(K,N_Z);
    Wcentroid_hull = Wcentroid;
    
    for k = 1:K
        c_k = C_new(:,k);
        %vert_temp = verts_solve(H,c_k);
        
        % Find Vertices of SSE B(W):
        vert_temp = lcon2vert(H,c_k,[],[],1e-8);
            VertSet{k} = vert_temp;
            
        % Q1: Take convhulln of vert_temp? 
        
        ind_temp = convhulln(vert_temp);
        Vert_chull_ind{k} = ind_temp;
        unique_v = unique(ind_temp);
        vert_chull_temp = vert_temp(unique_v,:);
        
        VertSetChull{k} = vert_chull_temp;
        
        % Q2: Then define directions and levels (G,c_in) ? for new set of
        % LP constraints instead of (H,c_new)? 
        
        [ G, cm, Geq, cmeq ] = vert2lcon(vert_chull_temp);
        
        % Normalized/direction vector that are normal to each face
        G_normal = [ G; Geq ];
        lengthG = sqrt(sum(G.^2, 2));
        
        if sum(lengthG) ~= numel(lengthG)
           lengthG = repmat(lengthG,1,N_Z);
           G_normal = G_normal ./ lengthG;
        end
        
        Gsearch{k} = G_normal;
        
        % Levels for each face
        c_normal{k} = [ cm; cmeq ];
        
        % Q3: If we go with Q2, then we can find CENTROID! 
        
        % Roughly, using center of containing hypersphere:
        [Wcentroid(k,:),~] = sphereFit(vert_chull_temp);
        
        % With some finesse, using the vertex points of convex hull to
        % calculate a weighted center of convex polytope:
        
         
         Wcentroid_hull(k,:) = centroid2( vert_temp(unique_v,:) );
    end
    

%%  Show and tell:
    figure('name','Original SSE Extreme Points $(H,c)$')
    
    for k = 1:K
       subplot(K/4,K/4,k)
       plot3(VertSet{k}(:,1), VertSet{k}(:,2), VertSet{k}(:,3), 'o', ...
                    Wcentroid_hull(k,1), Wcentroid_hull(k,2), ...
                                                Wcentroid_hull(k,3), 'or')
       campos([2 13 10])
       xlabel('w_1')
       ylabel('w_2')
       zlabel('w_3')
    end
    
%%  Show and tell:
    figure('name','Vertices of Convex-hull co(W)')
    
    for k = 1:K
       subplot(K/4,K/4,k)
       plot3(VertSetChull{k}(:,1), ...
                VertSetChull{k}(:,2), ...
                    VertSetChull{k}(:,3), 'o', ...
                    Wcentroid_hull(k,1), Wcentroid_hull(k,2), ...
                                            Wcentroid_hull(k,3), 'or') 
        campos([2 13 10])
        xlabel('w_1')
        ylabel('w_2')
        zlabel('w_3')
    end
    
    %% Interpolated data plot
    figure('name','Convex-hull co(W)')
    for k = 1:K
        x = VertSet{k}(:,1);
        y = VertSet{k}(:,2);
        z = VertSet{k}(:,3);
        ind_temp = Vert_chull_ind{k};
        
          
            subplot(K/4,K/4,k)
            hold on
                plot3(Wcentroid_hull(k,1), Wcentroid_hull(k,2), ...
                                Wcentroid_hull(k,3),... 
                                               'or', 'MarkerFaceColor','r') 
                h = trisurf(ind_temp,x,y,z, 'edgecolor',[0.5 0.5 0.5]);
                alpha(0.25)
                colormap(bone)
                campos([2 13 10])
                %camlight;
                %lighting gouraud;
            hold off
            box off
            set(h,'FaceLighting','gouraud',...
                        'FaceColor','interp',...
                            'AmbientStrength',0.9)
            light('Position',[1 0 0],'Style','infinite');
            set(gca,'color','none')
            xlabel('w_1')
            ylabel('w_2')
            zlabel('w_3')
    end

    
%% NEW SHIT ....

%% ACTION PROFILE SETS

    % Agent Action set
    a_j = linspace(0.01, amax(model), Nagrid)'; % j's action set
    %A_j = a_j;

    agrid = a_j;
    for i = 1 : N_Z-1
        agrid = gridmake(agrid,a_j);
    end

    ProfileSetA = agrid;
    NA = size(ProfileSetA,1);


    % Government action set
    b_u = linspace(0.01,model.MBAR,Nbgrid)';
    b_e = linspace(-model.WAGE*0.5, min(-0.1,model.MBAR),...
                                                     Nbgrid)';

    %b_u = linspace(0.0,self.MBAR,self.Nbgrid)';
    %b_e = linspace(-self.WAGE*0.25, min(-0.1,self.MBAR),self.Nbgrid)';
    bsu = b_u;
    for i = 1:N-1
        bsu = gridmake(bsu,b_u);
    end

    bse = b_e;

    for i = 1:M-1
        bse = gridmake(bse,b_e);
    end

    ProfileSetB = gridmake(bsu,bse); % requires GRIDMAKE: COMPECON
    NB = size(ProfileSetB,1);

    model.ProfileSetB = ProfileSetB;
    model.ProfileSetA = ProfileSetA;
    
    parallel = 'off';
    P_a = TransProbA(model,parallel);
    
    delstar = (1-DELTA);
    

%% CONSTRUCT SAMPLE SSE PATH (INDEXES TO (a,b) in A x B
alpha = cell(T,1);
beta = cell(T,1);
wvec = beta;
lams = beta;
                              
% Storage
    val_a = zeros(1,NA);
    wnext_a = cell(NA,1);
    
% Initial conditions
    k_now = InPartElement(lambda_now, dt, DtriVerts);
    
    %% UNDER DEVELOPMENT
    X = VertSetChull{k_now};    % Set of extreme points of W(k_now)
    
    
    %% Maximal Government payoff given (lambda_now, X)
    [ vG_max, index_max ]  = max(lambda_now * X');
    w_now = X(index_max,:)';
    
        % Check w_now is contained in W(k_now)
        tri = delaunayn(X);         % 3D-Delaunay triangulations
        tn = tsearchn(X, tri, w_now'); 
        IsInside = ~isnan(tn);
        
        if IsInside
            sprintf('\n Initial  w_now is inside W(k), k = %0.5g\n', k_now)
        else
            error('STRAT_INNERAY:IsInside',...
                   'Initial  w_now NOT inside W(k)')
        end
        
    
    % Best SSE payoff vector over Q(k_now):
    %     w_now = wmax(k_now,:)';
    %     
    %     v_now = vG_k_max(k_now);

    
%% MAIN LOOP: COMPUTING STRATEGIES SUPPORTING (lambda_now, w_now)   
     
for t = 1 : T
    
    v_now = lambda_now * w_now;
   
    %% Work out (approximately) feasible set of policies b, B(lambda)
    
    Feas_temp = lambda_now * ProfileSetB';
    B_temp = ProfileSetB( (abs(Feas_temp) <= tolB), :);
    NB = size(B_temp,1);
        
    % Inner loops:
    for b_idx = 1 : NB             % Loop of feasible policies b in B_temp
        
        % Variable storage:
        val_b = zeros(NB,1);
        aoptindex_b = cell(NB,1);
        wnext_b = cell(NB,1);
        
        parfor a_idx = 1 : NA       % loop over agents' actions a in A
            
            lb_local = lb;
            ub_local = ub;
            param_local = param;
            vartype_local = vartype;
            
            ProfileSetA_local = ProfileSetA;
            P_a_local = P_a;
            B_temp_local = B_temp;
            
            %C_new_local = C_new;
            c_normal_local = c_normal;
            G_normal_local = Gsearch;
            
            delstar_local = delstar;
            
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
            DevP = DELTA*(P_ia_not - P_ia);
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

                                                        
            % Get continuation game state:
            
            lambda_next = lambda_now * P_now;
            
            % Check which Q(k) contains lambda_next:        
            k_next = InPartElement(lambda_next, dt, DtriVerts);
            
            % Get appropriate partitition element's correspondence:
            %c_next = C_new_local(:,k_next);
            G_next = G_normal_local{k_next};
            c_next = c_normal_local{k_next};
            
            % Solve LP problem for maximizers w_next that enforces current
            % action-policy pair (a,b):
            
            constraint = [  G_next;                      %... (1)
                            DevP;                   %... (2)
                            lambda_now * DELTA*P_now   ];        %... (3)
                                        
            bconstraint1 = [ c_next;                              %...  (1)
                            DevA                      ];          %...  (2)
            bconstraint2 =  v_now - lambda_now * Util_now ;    %...  (3)
                                          
                bconstraint = [ bconstraint1; bconstraint2 ] ;
                
                
            objective = lambda_next; %ones(1,N_Z);
            
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
                    
                    if status == 2     
                        %status_vec(n) = status; 
                        %warning('STRAT:status',...
                        %            '[2] Feasible but not optimal');
                        val_a(a_idx) = fval;
                    elseif status == 5
                        %status_vec(n) = status;
                        val_a(a_idx) = fval;
                    else
                        val_a(a_idx) = -inf;
                        %warning('STRAT:status',...
                        %                      'Solution not found');
                    end          
        end % EndFor a
                
        % Given b, maximize over the a's:
        
        [ vmax_temp, ind_temp  ] = max(val_a);
        
        if isempty(ind_temp)
            error('STRAT:ind_temp','No feasible solution found')
        end
        
        val_b(b_idx) = vmax_temp; 
        aoptindex_b{b_idx} = ind_temp;
        
        wnext_ab_temp = wnext_a{ind_temp};
        wnext_b{b_idx} = wnext_ab_temp;
        
        fprintf('t = %i \t b_idx = %i \t \n',...
                                                t, b_idx);
        %disp(wnext_ab_temp);
        
    end % EndFor b
    
    % Maximize over the feasible b's:
    
    [~, boptindex ] = max(val_b);
    
    if isempty(ind_temp)
            error('STRAT:boptindex','No feasible solution found')
    end
        
    w_next = wnext_b{boptindex};
        
    % Store profiles (a*, b*) at stage/date t:    
    alpha{t} = ProfileSetA(aoptindex_b{boptindex},:); % Index to row vector in A
    beta{t} = B_temp(boptindex,:);        % Index to row vector in B_temp
   
    % Update state and continuation payoffs
    lambda_next = lambda_now * P_a(:,:,aoptindex_b{boptindex});
    
    lambda_now = lambda_next;
    w_now = w_next;
    
    wvec{t} = w_now;
    lams{t} = lambda_now;
                           
end % EndFor t

% % Pack it all in:
strat.a = alpha;
strat.b = beta;
strat.w = wvec;
strat.l = lams;
                           
save('_mat/strat.mat', 'strat')

