% % SSE_COMPARE.M
% %
% % Calculate:
% % (1) CASE 1: Best and worst values sustainable in an SSE;
% % (2) CASE 2: Best steady-state commitment policy payoff;
% % (3) CASE 3: Best dynamic Ramsey commitment policy payoff, 
% % for every initial aggregate state (equivalently, game state).
% %     
% %
% % =======================================================================
% % (c) 2011-- T.Kam and R.Stauber.
% %
% % DEPENDENCIES:
% %
% %     * GNU Linear Programming Kit (ANSI C to MATLAB MEX)
% %     * CompEcon Toolbox (Miranda-Fackler)
% %         
% % Email: mortheus__at__gmail.com
% % =======================================================================
% % $Revision: 4.0.3 $  $Date: 2011/06/01 11:40:20 $ 
% %
% % See also AGENT, MAIN, GETINDEX, GLPK
close all
clear all
clc

UNIX = 'on';                    % UNIX users: 'on' ||  MS-Windows: 'off'

% Add path to directory:
    if strcmp(UNIX,'on') == 1   
        addpath(strcat(pwd,'/glpkmex/'));
        addpath(strcat(pwd,'/_func/'));
        addpath(strcat(pwd,'/_mat/'));
        addpath(strcat(pwd,'/_figures/'));
        addpath(strcat(pwd,'/eq_sphere_partitions/'));
    else
        addpath(strcat(pwd,'\glpkmex\'));
        addpath(strcat(pwd,'\_func\'));
        addpath(strcat(pwd,'\_mat\'));
        addpath(strcat(pwd,'\_figures\'));
        addpath(strcat(pwd,'\eq_sphere_partitions\'));
    end


% CASE 1: 

    % Load results from MAIN.M (SSE Outer Approximation):
    load('_mat/sse_iter438')
    
    vG_sse = sse_extreme(C_old, np, V_MAX, V_MIN);

% CASE 2: Optimal steady state with fixed policy commitment outcome

    % Load results from STEADYSTATE.M (Fixed policy DP problem):
    load ramsey_lr  
    k_ramsey = getindex(lambda_ramsey_ss_max, D);


% CASE 3: Optimal Ramsey commitment outcome
    clear C_old
    load('_mat/ramsey_iter239');
    vG_ramsey = ramsey_extreme(C_old, np, V_MAX, V_MIN);


% FIGURES
% -------------------------------------------------------------------------
figure('name','SSE worst/best social values at each state k')

subplot(1,2,1)
    hold on

    plot(   1:K, vG_sse(:,1), 'or', ...
            1:K, vG_sse(:,2), 'ob', ...
            1:K, vG_ramsey(:,2), 'xm', ...
            k_ramsey, vrmax, 'rs' )

    
    legend( 'SSE v_{G}^{min}(\lambda_k)', ...
            'SSE v_{G}^{max}(\lambda_k)', ...
            'Ramsey v_{G}^{max}(\lambda_k)', ...
            'L/R Commitment, v_{G}^{comm}' )
                               
    text(k_ramsey, vrmax, sprintf('(%i,%0.2g)',k_ramsey, vrmax),...
                'HorizontalAlignment','right')
    
    legend boxoff
    
    hold off                               
    xlabel('State Index, k')
    ylabel('Current social utility')
    
    subplot(1,2,2)
        stem([-1:-1:-N,1:M],lambda_ramsey_ss_max,'MarkerFaceColor','r')
        xlim([-2 2])
        xlabel('Individual State')
        ylabel(['\lambda (k), k = ', int2str(k_ramsey)])

    hgsave('_figures/sse_compare_fig1')
    print -depsc _figures/sse_compare_fig1

        
% FIGURE 2:
figure('name','SSE best social values at each state k')

    subplot(1,2,1)
        hold on

        plot(   k_ramsey, vG_sse(k_ramsey,2), 'ob', ...
                k_ramsey, vG_ramsey(k_ramsey,2), 'xm', ...
                k_ramsey, vrmax, 'rs' )


        legend( 'SSE v_{G}^{max}(\lambda_k)', ...
                'Ramsey v_{G}^{max}(\lambda_k)', ...
                'L/R Commitment, v_{G}^{comm}' )

        text(k_ramsey, vrmax, sprintf('(%i,%0.2g)',k_ramsey, vrmax),...
                    'HorizontalAlignment','right')

        %legend boxoff

        hold off                               
        xlabel('State Index, k')
        ylabel('Current social utility')

    subplot(1,2,2)
        stem([-1:-1:-N,1:M],lambda_ramsey_ss_max,'MarkerFaceColor','r')
        xlim([-2 2])
        xlabel('Individual State')
        ylabel(['\lambda (k), k = ', int2str(k_ramsey)])
 
    hgsave('_figures/sse_compare_fig2')
    print -depsc _figures/sse_compare_fig2

MSZ = 10; % Scale up MarkerSize from default 6pt to MSZpt.

% FIGURE 3: Comparing best SSE and Ramsey over some subset of state space
j = 2;  j = max(j,1); % lowest k. Note min(k) = 1
J = 15; J = min(J,K); % largest k. Note: max(k) = K;

 vcomp = [vG_sse(:,2), vG_ramsey(:,2)];
 
figure('name','Comparing best SSE and Ramsey')
    grid on
    plot(j:J, vcomp(j:J,1),'ob', j:J, vcomp(j:J,2),'xm','MarkerSize', MSZ)
    xlabel('State Index, k')
    ylabel('Current social utility')
    
    hgsave('_figures/sse_compare_fig3')
    print -depsc _figures/sse_compare_fig3
    
% FIGURE 4: Comparing best SSE and Ramsey over some subset of state space
j = 16;  j = max(j,1); % lowest k. Note min(k) = 1
J = 17; J = min(J,K); % largest k. Note: max(k) = K;


figure('name','Comparing best SSE and Ramsey')
    grid on
    plot(j:J, vcomp(j:J,1),'ob', ...
            j:J, vcomp(j:J,2),'xm', 'MarkerSize', MSZ)
    xlabel('State Index, k')
    ylabel('Current social utility')
    
    hgsave('_figures/sse_compare_fig4')
    print -depsc _figures/sse_compare_fig4  
    
% FIGURE 4: Comparing best SSE and Ramsey over some subset of state space
j = 18;  j = max(j,1); % lowest k. Note min(k) = 1
J = 19; J = min(J,K); % largest k. Note: max(k) = K;

figure('name','Comparing best SSE and Ramsey')
    grid on
    plot(j:J, vcomp(j:J,1),'ob', j:J, vcomp(j:J,2),'xm', 'MarkerSize', MSZ)
    xlabel('State Index, k')
    ylabel('Current social utility')
    
    hgsave('_figures/sse_compare_fig5')
    print -depsc _figures/sse_compare_fig5