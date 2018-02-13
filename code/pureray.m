function [ l_star, aopt_idx, bopt_idx, L_m ] = pureray(model, theta_ray,...
                            lambda_now, k_now, B, P_a, ...
                            wko, dt, DtriVerts, G_normal, c_normal, pival)

    %% PURERAY.M
    %
    % Given current state lambda_now in partition element k_now, and,
    % for a given theta_ray direction, compute (if it exists) a scalar
    % length l(a,b) where the (a,b) profile could support the initial
    % value centroid(k_now) + theta_ray*l(a,b) as an SSE. The output L_m
    % records these results l(a,b) for every profile (a,b). Note some or
    % all l(a,b) could be -inf valued, meaning that there is no
    % profile-and-extreme-point tuple,(a,b,theta_ray), that sustains a
    % pure-strategy SSE in the *particular direction* theta_ray.
    %
    % Inputs:
    % -------
    % model: Model objects
    % theta_ray: Current ray direction (candidate extreme point)
    % lambda_now: 1 x NZ row vector of agent distribution (current)
    % k_now: scalar index of partition element Q_k containing lambda_now
    % B: NB x NZ array of budget-balancing policies at state lambda_now
    % P_a: Array of transition matrices (each indexed by a profile "a")
    % wko: NZ x 1 array---centroid of W_in(k)
    % dt, Dtriverts: Delaunay triangularization of D
    % G_normal, c_normal: directions and levels of W_in facets
    % pival: K x 1 array of punishment values at each k = 1,..,K
    %
    % Outputs:
    % --------
    % L_m: NA x NB array of lengths l(a,b)
    % l_star: = max(max(L_m))
    % aopt_idx: = argmax(L_m,[],1) (index position)
    % bopt_idx: = argmax(max(L_m)) (index position)
    %
    % (c) 2012, 2013, Timothy Kam. Email: tcy.kam__at__gmail.com
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Use subject to GNU LGPL licensing terms. Cite this header and
    % author completely in subsequent re-use and modifications.
    % =====================================================================
    % $Revision: 5.0.0 $  $Date: 2015/12/31 08:45:20 $

    A = model.ProfileSetA;
    NA = size(A,1);
    NB = size(B,1);
    N = model.N;
    delta = model.DELTA;
    deltstar = model.DELSTAR;
    N_Z = model.NS;

    % hack
    %H = model.H;
    %c = model.c;

    % Centroid of W_in(k):
    %     wko = centroids(k_now,:)';

    % Punishment value at k_now:
    pi_k = pival(k_now);

    L_m = -inf(NA,NB); % Preallocation with -inf
    for b_idx = 1:NB
        for a_idx = 1:NA
            % State Transition for a at a_idx:
                lambda_next = lambda_now * P_a(:,:,a_idx);
                k_next = InPartElement(lambda_next, dt, DtriVerts);
            % Get appropriate partitition element's correspondence:
                G_next = G_normal{k_next};
                c_next = c_normal{k_next};

                % hack
                %G_next = H;
                %c_next = c(:,k_next);

            % Current payoff from a: effort disutility
                Util_a = -phia(model, A(a_idx,:));
                Util_ia = repmat(Util_a, NA-1, 1);
            % Current payoff from b:
                b_current_feasible = B(b_idx,:);
                cons = b_current_feasible;
                cons(N+1:end) = cons(N+1:end) + model.WAGE;
                Util_b = ucons(model, cons);
            % Total current payoff: (column vector!)
                Util_now = deltstar*( Util_a + Util_b )';

            % Constraint (1): Consistency w.r.t. W_in(k_next)
            % ------------------------------------------------------
            Consistency_lhs = G_next*(1/delta)*theta_ray;
            Consistency_rhs = c_next-G_next*(1/delta)*(wko-Util_now);

            % Constraint (2): Government's IC
            % ------------------------------------------------------
            GovIC_lhs = - lambda_now*theta_ray;
            GovIC_rhs = -pi_k + lambda_now*wko;

            % Constraint (3): Agents' CE conditions (IC)
            % ------------------------------------------------------
                % All other a# \neq a, and induced P(a#) matrices:
                Index_A = 1:NA;                     % indexes to { a }
                ia_Not = setdiff(Index_A, a_idx);
                                                    % All a# \neq a
                P_ia_Not = P_a(:,:,ia_Not);   % All P(a#) \neq P(a),
                                                    % N x N x card(A)-1
                P_ia_not = reshape(P_ia_Not,N_Z*(NA-1),N_Z);
                                                % Tile this too!

                % Current payoffs from all deviations from a:
                Util_ia_Not = -phia(model, A(ia_Not,:) );

                % Transitions
                P_now = P_a(:,:,a_idx);
                P_ia = repmat(P_now,NA-1,1);

                % Agents' incentive constraints at current {k,l,a}
                DevP = delta*(P_ia_not - P_ia);
                DevA = deltstar*(Util_ia - Util_ia_Not);
                [nr, nc ] = size(DevA);
                DevA = reshape(DevA', nr*nc, 1);
            CE_lhs = DevP*(1/delta)*theta_ray;
            CE_rhs = DevA - DevP*(1/delta)*(wko - Util_now);

            % ------------------------------------------------------------
            % Solve LP problem for maximizers w_next that enforces current
            % action-policy pair (a,b).
            % Linear objective in length, l:
            objective = 1;
            % NOTE: (1)-(3): weak inequality constraints:
            Constraint_lhs = [  Consistency_lhs;                   %... (1)
                                GovIC_lhs;                         %... (2)
                                CE_lhs           ];                %... (3)

            Constraint_rhs = [  Consistency_rhs;                   %... (1)
                                GovIC_rhs;                         %... (2)
                                CE_rhs           ];                %... (3)

            % Constraint types: weak inequalities
            lb_local = 0; % l > 0
            ub_local = [];
            ctype_local = repmat('U', size(Constraint_rhs,1),1);
            vartype_local = 'C';

            %% GLPK Settings:
            sense = -1;     % Max (-1) or Min (+1) problem
            param.msglev=0; % Output GLPK messages
            param.save=0;   % Set save options

            [lopt,~,status,~] = glpk( objective, Constraint_lhs, ...
                                       Constraint_rhs,lb_local,ub_local,...
                                       ctype_local,vartype_local,...
                                       sense,param);
            if status == 5
                L_m(a_idx,b_idx) = lopt;
            end
        end
    end
    [lb_temp, aopt_idx] = max(L_m,[],1); % Max over a_idx (rows)
    %aopt_idx
    [l_star, bopt_idx] = max(lb_temp); % Max over b_idx
    aopt_idx = aopt_idx(bopt_idx);
end
