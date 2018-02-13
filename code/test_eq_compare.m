% test_eqcompare.m

clear all
close all
clc


%% COMPUTE THREE POLICY CASES:

    %% CASE 1: SSE payoff correspondence:
    load('_mat/sse2014_iter162.mat')

        L = size(lambda_new, 3);
        K = size(lambda_new, 4);

        % lambda_new and w_new are (1 x 3 x L x K ) arrays
        vG_k = sum( lambda_new .* w_new, 2 );
        vG_k = reshape(vG_k, L, K);
        vG_k_max = max(vG_k,[],1);
        vG_k_min = min(vG_k,[],1);


    %% CASE 2: Ramsey payoff correspondence:
     load('_mat/ramsey2014_iter48.mat')

        % lambda_new and w_new are (1 x 3 x L x K ) arrays
        vR_k = sum( lambda_new .* w_new, 2 );
        vR_k = reshape(vR_k, L, K);
        vR_k_max = max(vR_k,[],1);
        %vR_k_min = min(vR_k,[],1);
        
    
    %% CASE 3: Steady-state optimal policy:

    load('_mat/ramsey_lr01.mat')
    
    % Check where lambda_ramsey_ss_max belongs to in simplex partition:
    
    kss = InPartElement(lambda_ramsey_ss_max, dt, DtriVerts);
    
%     cycle = [ 1 2 3 1 ];
%     
%     x = lambda_ramsey_ss_max(1);
%     y = lambda_ramsey_ss_max(2);
%     
%     in = zeros(K,1);
%     on = in;
%     
%     for k = 1:K
%         Q_k = DtriVerts(dt.Triangulation(k,:),:);
%         Q_k = Q_k(cycle,:);
%         xv = Q_k(:,1);
%         yv = Q_k(:,2);
%         
%         [in(k), on(k)] = inpolygon(x,y,xv,yv);
%         
%     end
%     
%     % Find index to partition element Q(k) containing lambda_ramsey_ss_max:
%     loc1 = find(in == 1);
%     loc2 = find(on == 1);
%         
%     if ~isempty(loc1) && isempty(loc2)
%         kss = loc1;   
%     elseif isempty(loc1) && ~isempty(loc2)
%         kss = loc2;
%     elseif ~isempty(loc1) && ~isempty(loc2)
%         u = rand(1);
%         if u < 0.5
%             kss = loc1;
%         else
%             kss = loc2;
%         end
%     else
%         warning('Empty kss!')
%     end
  
        
        


%% PLOT RESULTS:
    figure('name', 'SSE v(G) bounds')
        hold on
        plot(1:K,vG_k_max,'sk', 1:K, vG_k_min,'*k', ... % SSE payoffs
                1:K, vR_k_max,'db')                      % Ramsey
        plot(1:K, vrmax*ones(size(1:K)), 'o-r') % Steady State 
        hold off
        
        legend( 'Best SSE', ...
                'Worst SSE', ...
                'Ramsey', ...
                'Optimal Steady State', ...
                'Location','Best' );
        xlabel('Index k of Partition Element Q(k)')
        ylabel('Average total welfare')
        print -depsc _figures/payoff_compare
