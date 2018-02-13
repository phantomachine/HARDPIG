% % STEADYSTATE.M
% %
% % USAGE: Main script for solving agents' best responses to all policy  
% % vectors (fixed). Then choose "best" fixed policy vector that maximizes
% % steady-state welfare.
% %
% % =======================================================================
% %     (c) 2011-- T.Kam and R.Stauber.
% %
% % DEPENDENCIES:
% %
% %     * KSMOD2 class for model and SPE solution.
% %     * GNU Linear Programming Kit (ANSI C to MATLAB MEX)
% %     * CompEcon Toolbox (Miranda-Fackler)
% %         
% % Email: mortheus__at__gmail.com
% % =======================================================================
% % $Revision: 4.0.3 $  $Date: 2014/03/08 12:38:20 $  
% %
% % See also SSE_EXTREME, RAMSEY, KSMOD02

clear all
close all
clc

LOAD_OLD = 0;

SUBPLOT = 'off';

%% ------------------------------------------------------------------------
%% STEP -1: MACHINE SPECIFIC SETTINGS
%% ------------------------------------------------------------------------


% Parallel Computing option (else, set N_PROCESSOR = 0)
N_PROCESSOR = 2;


%% ------------------------------------------------------------------------
%% STEP 0: MODEL PARAMETERIZATION
%% ------------------------------------------------------------------------

SIGMA = 1;                % CRRA u(.) parameter (linear = 0; log = 1)
PHI = 0.0;                  % Disutility of effort parameter (linear = 0)
M = 1;                    % "Max duration" of employment, M
N = 2;                    % "Max duration" of unemployment, N
WAGE = 1;                 % fixed common wage for employed
REP_RATIO = 0.8;          % "Maximal replacement ratio"
DELTA = 0.9615;           % common discount factor: Annual rate = 4%
NGRID_P = 8;              % # grid points on [0,1] in SSE code, redundant here

% NOTES: Dependant parameters set in KSMOD2.M
% -------------------------------------------
% MBAR = WAGE*REP_RATIO;    % exogenous upperbound on benefit payout
% DELSTAR = (1 - DELTA);    % For averaging payoffs

%% ------------------------------------------------------------------------
%% STEP 1: CREATE INSTANCE OF CLASS: KSMOD02
%% ------------------------------------------------------------------------

model = ksmod03ss(SIGMA,PHI,M,N,WAGE,REP_RATIO,DELTA); % Given STEP 0
%ksmod03(SIGMA,PHI,M,N,WAGE,REP_RATIO,DELTA)
%% ------------------------------------------------------------------------
%% STEP 2: STORE AGENT'S BEST ACTION FUNCTION OUTSIDE
%% ------------------------------------------------------------------------

N_BGRID = size(model.ProfileSetB,1);
N_Z = model.NS;

% Initialize memory space for (v,a):
v = zeros(N_BGRID, N_Z);
lambda_ss = v;
a = v;
surplus = zeros(N_BGRID,1);

display3 = model.Display3;
b_set = model.ProfileSetB;

% Use Parallel Computing Toolbox and Howard's Improvement Algorithm:
%matlabpool(N_PROCESSOR);
if LOAD_OLD ~= 1
    
    if matlabpool('size') > 0
    matlabpool close force
    end
    
    matlabpool(N_PROCESSOR);

    tic

    fprintf('\nPLEASE WAIT: Solving for (v,a) - iteration over b:\n\n');
    parfor count = 1 : N_BGRID

        if  display3 == 1
           fprintf('No:%i of %i\n', count, N_BGRID);
        end

        % Current policy function of personal states: b(-N,...,-1,1,...,M):
        b = b_set(count,:);

        % f(-1,...,-N,1,...,M) where f := [v_opt,a_opt]:
        [v_opt,a_opt] = howard(model,b);

        % Store value and best reply functions as functions of b's:
        v(count,:) = v_opt ; % Store counter NB in column #1
        a(count,:) = a_opt ; % for sorting later on if needed.

        au = a_opt(1:N);
        ae = a_opt(N+1:end);
        P = markovmat(model, au, ae);

        % Solve for ergodic distribution
        z = zeros(1,N_Z);
        z(end) = 1;
        PMI = P - eye(N_Z);
        PMI(:,end) = ones(N_Z,1);
        lambda = z / PMI;

        lambda_ss(count,:) = lambda;

        surplus(count) = lambda * b';
    end

    toc

    runTime = toc;
    avgTime = toc/N_BGRID;
    
    matlabpool close

    fprintf('DONE ... Solving for (v,a)\n\n');

    fprintf('START ... Finding feasible b and induced (v,a)\n\n');

    % Get logical index set of SURPLUS that satisfies small nonegativity
    feas_condition = (surplus >= -0.15 & surplus <= 0);

    % Map into payoff vectors
    v_ss = v(feas_condition,:);

    % Map into transfers
    b_ss = b_set(feas_condition,:);
    
    % Map into actions
    
    a_ss = a(feas_condition,:);

    % Map into distribution
    
    lambda_ramsey_ss = lambda_ss(feas_condition,:);

    % Given v_ss and Lambda_ramsey_ss sets, calculate Ramsey payoffs
    v_ramsey = sum( lambda_ramsey_ss .* v_ss, 2);

    % Now choose the Ramsey payoff consistent with (near)-zero surplus
    feas_surplus = max( surplus(feas_condition) );

    [vrmax,k] = max(v_ramsey (surplus(feas_condition) == feas_surplus) );

    % Best Ramsey (approximate) payoff
    vss_max = v_ss(k,:);
    
    % Best Ramsey (approximate) transfer
    bss_max = b_ss(k,:);
    
    % Best Ramsey private actions
    ass_max = a_ss(k,:);
    
    % Best Ramsey payoff
    lambda_ramsey_ss_max = lambda_ramsey_ss(k,:);
    

    % SAVE RESULTS FOR SSE_EXTREME.M
    save('_mat/ramsey_lr01.mat');

elseif LOAD_OLD == 1
    load('_mat/ramsey_lr01.mat');
end

% DISPLAY RESULTS
fprintf('\nMax. Ramsey steady state social value = %g\n\n',vrmax)
disp('Induced steady state distribution:')

disp(lambda_ramsey_ss_max)

if strcmp(SUBPLOT,'on')
    figure('name','Steady-state Ramsey Outcomes')

        axis tight

        subplot(2,2,1)
        stem([-1:-1:-N,1:M],lambda_ramsey_ss_max,'MarkerFaceColor','r')
        xlim([-2 2])
        grid on
        xlabel('Z')
        ylabel('Probability')
v_
        subplot(2,2,2)
        stem([-1:-1:-N,1:M],bss_max,'MarkerFaceColor','r')
        grid on
        xlabel('Z')
        ylabel('Tranfers')

        subplot(2,2,3)
        stem([-1:-1:-N,1:M],ass_max,'MarkerFaceColor','r')
        grid on
        xlabel('Z')
        ylabel('Effort')

        subplot(2,2,4)
        stem([-1:-1:-N,1:M],vss_max,'MarkerFaceColor','r')
        grid on
        xlabel('Z')
        ylabel('Indirect Utilities')

        print -depsc _figures/fixed_ss1
else
    
    % Reorder the -j outcomes:
    Zdomain = [-N:-1,1:M];
    lambda = zeros(1,M+N);
        lambda(1:N) = lambda_ramsey_ss_max(1:N);
        %lambda(1:N) = lambda_ramsey_ss_max(N:-1:1);
        lambda(N+1:end) = lambda_ramsey_ss_max(N+1:end);
        
    wss = zeros(1,M+N);
        %wss(1:N) = vss_max(N:-1:1);
        wss(1:N) = vss_max(N:-1:1);
        wss(N+1:end) = vss_max(N+1:end);
        
    bss = zeros(1,M+N);
        %bss(1:N) = bss_max(N:-1:1);
        bss(1:N) = bss_max(1:N);
        bss(N+1:end) = bss_max(N+1:end);
        
    ass = zeros(1,M+N);
        %ass(1:N) = ass_max(N:-1:1);
        ass(1:N) = ass_max(1:N);
        ass(N+1:end) = ass_max(N+1:end);
        
        
    figname = { 'Steady-state Optimal Distribution', ...
                'Steady-state Optimal a Profile', ...
                'Steady-state Optimal b Profile', ...
                'Steady-state Optimal w Profile'         };
            
    Ylabel_data = {     'Proportion', ...
                        'Action a(j)', ...
                        'Policy b(j)', ...
                        'Payoff w(j)',                  };
            
    filename_data = {   'fixed_ss_dist', ...
                        'fixed_ss_a', ...
                        'fixed_ss_b', ...
                        'fixes_ss_w'                    };
                    
                    
    data = [ lambda; ass; bss; wss ];
    for n = 1:numel(figname)
        
        figure('name',figname{n})
            bar(data(n,:), 0.15,'c','EdgeColor','b');
            for i = 1:N+M
               if abs(data(n,i)) > 0.1
                   heightscale = 0.8;
               elseif abs(data(n,i)) <= 0.1
                   heightscale = 1.5;
               end

               text(i,data(n,i)*heightscale,...
                           sprintf('%3.3g', data(n,i)),'Color','k' );
            end
            
            set(gca,'XTickLabel',{Zdomain})
            grid off
            xlabel('Z')
            ylabel(Ylabel_data{n})
            
            currfilename = strcat('_figures/',filename_data{n});
            print('-depsc', currfilename);
            hgsave(currfilename);
    end      
            

    
end