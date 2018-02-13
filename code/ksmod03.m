% % Class -- KSMOD03.M
% % =======================================================================
% %     (c) 2013-- T.Kam
% %         Defines class for model and SPE solution.
% %
% %         Email: tcy.kam@gmail.com
% % =======================================================================
% % $Revision: 4.1.0 $  $Date: 2013/09/08 11:18:00 $

classdef ksmod03 < handle

    properties
        SIGMA;      % CRRA parameter (linear = 0; log = 1)
        PHI;        % Disutility of effort parameter (linear = 0)
        M;          % "Max duration" of employment, M
        N;          % "Max duration" of unemployment, N
        WAGE;       % fixed common wage for employed
        REP_RATIO;  % Max. Replacement ratio
        DELTA;      % Common discount factor

        MBAR;       % upper bound on benefit payout
        DELSTAR;    % 1 - DELTA
        SBAR;       % upper bound on surplus/saving

        Z;          % States of unemployment and employment duration
                    % Note: State space: R x D(Z)
        NS;         % dim(D(Z)) - no. of relevant individual states
        NP;         % dim(payoff) (workers states + government)

        A_j;        % Discretized space for a (agent action profile space)
        ProfileSetA;% Discretized space for A_j^NS (action profile space)
        ProfileSetB;% Grid array for all b's (govt action profile space)

        %NGRID_P;    % #elements in Grid approximating [0,1]
        Nagrid = 2; % No. of grid points for AGRID
        Nbgrid = 2; % No. grid points for each b(i), i = -N,...,-1,1,...,M
        Results;    % Object to store results

        NA;         % # private agent action profiles, card(A)
        NB;         % # government action profiles, card(B)
    end

    properties (Access = public)
        H;
        c;
    end

    properties (Constant)
        TOL1 = 0.001;
        TOL2 = 1e-8;
        MAXITER1 = 500;
        MAXITER3 = 500;
        MAXITER2 = 500;
        Options = optimset('display','on');
        Display1 = 0;
        Display2 = 0;
        Display3 = 1;
    end

    methods
        function self = ksmod03(SIGMA,PHI,M,N,WAGE,REP_RATIO,DELTA)
            % KSMOD02:
            % Set up various model parameters and state spaces. Uses
            % GRIDMAKE from COMPECON Toolbox by Miranda and Fackler
            %
            % See also GRIDMAKE

            % Parameters and bounds
            self.SIGMA = SIGMA;
            self.PHI = PHI;
            self.M = M;
            self.N = N;
            self.WAGE = WAGE;
            self.REP_RATIO = REP_RATIO;
            self.DELTA = DELTA;
            self.MBAR = WAGE*REP_RATIO;
            self.DELSTAR = (1 - DELTA);
            self.SBAR = WAGE/(1-DELTA);

            % # Probability realizations
            % self.NGRID_P = NGRID_P;

            % Individual state space
            self.Z = [-(1:N),(1:M)];
            self.NS = length(self.Z);
            self.NP = self.NS + 1;


            % Agent Action set
            a_j = linspace(0.01, amax(self), self.Nagrid)'; % j's action set
            self.A_j = a_j;

            agrid = a_j;
            for i = 1 : self.NS-1
                agrid = gridmake(agrid,a_j);
            end

            self.ProfileSetA = agrid;
            self.NA = size(self.ProfileSetA,1);


            % Government action set
            b_u = linspace(0.01,self.MBAR,self.Nbgrid)';
            b_e = linspace(-self.WAGE*0.5, min(-0.1,self.MBAR),...
                                                             self.Nbgrid)';

            bsu = b_u;
            for i = 1:N-1
                bsu = gridmake(bsu,b_u);
            end

            bse = b_e;

            for i = 1:M-1
                bse = gridmake(bse,b_e);
            end

            self.ProfileSetB = gridmake(bsu,bse); % requires GRIDMAKE: COMPECON
            self.NB = size(self.ProfileSetB,1);

        end

        function IndexProfileSet = ProfileSetIndex(self)
            % ProfileSetIndex.M
            % Indexing array to action profiles (a, b)

            %NA = size(self.ProfileSetA,1);
            IndexA = (1:self.NA)';

            %NB = size(self.ProfileSetB,1);
            IndexB = (1:self.NB)';

            IndexProfileSet = gridmake(IndexA, IndexB);
        end

        function uc = ucons(self,cons)
            % UCONS.M
            % Two cases for the CRRA class of utility function.

            sigma = self.SIGMA;

            if sigma == 1
                uc = log(cons + 1e-12);
            elseif sigma < 0
                warning('UCONS:Error','SIGMA must be nonnegative')
                return
            else
                uc = (cons + 1e-12).^(1-sigma)/(1-sigma);
            end
        end

        function fi = phia(self,a)
            % PHIA.M
            % Two cases for the CRRA class of utility function.

            phi = self.PHI;

            %if phi == 1
            %    fi = log(a + 1e-12);
            if phi < 0
                warning('PHIA:Error','PHI must be nonnegative')
                return
            else
                fi = a.^(1+phi)/(1+phi);
            end
        end

        function a_up = amax(self)
            % AMAX.M:
            % upper bound on payoffs: a = 0, cons = WAGE + MBAR
            V_MAX = ucons(self, (self.WAGE + self.MBAR)); %/(1-self.DELTA);

            % solve for upper bound on effort, astar
            % require Symbolic Math Toolbox
            syms x
            astar = solve(self.DELTA*exp(-x)*V_MAX ...
                                          - (1 - self.DELTA)*(x^self.PHI));
            a_up = double(astar);
        end

        function prob_e = pe(self,a)
            % PE.M:
            %
            % Input
            %   a           :   NGRID_A x M matrix of positive reals
            %   j           :   1 x M vector of +ve integers
            % Output
            %   prob_u      :   Transition probabilities at each (a,j)

                m = self.M;
                prob_e = 1 - 1./(repmat((1:m),size(a,1),1).*exp(a));
        end

        function prob_u = pu(self,a)
            % PU.M:
            %
            % Input
            %   a           :   NGRID_A x N matrix of positive reals
            %   j           :   1 x N vector of -ve integers
            % Output
            %   prob_u      :   Transition probabilities at each (a,j)
                n = self.N;
                prob_u = (1 - 1./exp(a))./repmat(-(-n:-1),size(a,1),1);
        end

        function vnext = conval(self, p_ue, p_ee, v)
            % CONVAL.M
            % Specific continuation value.
            %
            % INPUT:
            %      self    :    class
            %      p_ue    :    NA x N matrix of t.p.f. for -j's
            %      p_ee    :    NA x N matrix of t.p.f. for +j's
            %      v       :    NA x (N+M) matrix, value function guess
            %
            % OUTPUT:
            %      vnext   :    NA x 1 vector, continuation values
            %
            % See also HOWARD, PE, PU

            n = self.N;

            np = size(p_ue,1);

            % Unemployed continuation value
            conval_u = p_ue*v(n+1) ...
                           + (1 - p_ue).*repmat(v([1, 1:n-1]),np,1);

            % employed's expected continuation value:
            vup = [v(n+2:end), v(end)];     % continuation value +1 step up
            conval_e = p_ee.*repmat(vup,np,1) + (1-p_ee).*v(n);

            vnext = [conval_u, conval_e];
        end

        function [v, a] = howard(self,taxtrans)

            % Unpack parameters
            delta = self.DELTA;
            delstar = self.DELSTAR;
            wage = self.WAGE;

            agrid = self.A_j;
            nstate = self.NS;
            n = self.N;
            m = self.M;
            na = self.NA;

            % Consumption set
            c = [zeros(1,n), wage*ones(1,m)] ...
                                        + taxtrans; % Budget constraint
            uc1 = ucons(self,c);

            uc = repmat(uc1,na,1);
            %uc = ucons(self,c);

            % Initialize:

            v = zeros(1,nstate);        % Value function v(i), i = 1,...,ns
            a = self.AGRID(2)*ones(1,nstate);

            amat = repmat(agrid,1,nstate);       % decision space

            p_ue1 = pu(self, a(:,1:n));
            p_ee1 = pe(self, a(:,n+1:n+m));

            p_ue = pu(self, amat(:,1:n));   % t.p.f. u to e
            p_ee = pe(self, amat(:,n+1:n+m)); % t.p.f. e to e+


            iter1 = 0;
            iter2 = 0;

            crit1 = 1;
            crit2 = 1;

            % Howard's Policy improvement
            while crit1 > self.TOL2 && iter1 < self.MAXITER3
                iter1 = iter1 + 1;

                while crit2 > self.TOL2 && iter2 < self.MAXITER2

                    iter2 = iter2 + 1;

                    Tv_afix = delstar*(uc1 - phia(self,a)) ...
                                       + delta*conval(self,p_ue1,p_ee1,v);

                    crit2 = max(abs(Tv_afix - v));

                    if self.Display2 == 1
                        disp([iter2, crit2])
                    end

                    v = Tv_afix;
                end

                % One-shot deviation: Two-period max problem, update policy

                [Tv, opt_index] = max(delstar*(uc - phia(self,amat)) ...
                                    + delta*conval(self,p_ue,p_ee,v),[],1);
                                                   % maximize along rows: a


                crit1 = max(abs(Tv - v));

                v = Tv;
                a = amat(opt_index);

                if self.Display1 == 1
                    disp([iter1, crit1])
                end
            end
        end

        function P = markovmat(self,au,ae)
            % MARKOVMAT.M
            % This function defines the Markov matrix for the dynamic
            % game with N+M finite states. Note: -N is minimal negative
            % integer state and +M is maximal positive integer state.
            %
            % EXAMPLE:
            % In our particular setup there will be a set of states
            % representing unemployment states (duration):
            %
            % S1 = { -N, ..., -1 }
            %
            % Each action vector (au, ae)
            %
            % and a set of states representing employment states:
            %
            % S2 = { 1, 2, ..., M }
            %
            % The markov matrix P will have the first N rows for S1 and the
            % remaining M rows will be for S2.
            %
            % Input:
            % au        1 x N vector of actions each for each state in S1
            % ae        1 x M vector of actions each for each state in S2
            %
            % Output:
            % P         ns x ns Markov matrix for current actions; ns = N+M
            %
            % See also PE, PU

            m = self.M;
            n = self.N;

            % Initialize Markov matrix component blocks
            P12 = zeros(n,m);
            P21 = zeros(m,n);

            % Current states in S1: transition probs
            pUU = 1 - pu(self,au);
            pUE = 1 - pUU;

            % Current states in S2: transition probs
            pEE = pe(self,ae);
            pEU = 1 - pEE;

            % Fill in the blanks
            P11 = diag(pUU(1:end-1),1);
            P11(end,end) = pUU(end);
            P12(:,1) = pUE';
            P21(:,1) = pEU';
            P22 = diag(pEE(1:end-1),1);
            P22(end,end) = pEE(end);

            % Concantenate to make markov matrix P = P(au,ae)
            P = [ P11, P12;
                  P21, P22  ];

            % Re-sort the PUU trasitions to order (-N,..,-1)
            P(:,1:n) = P(:,n:-1:1);
        end

        function dcoeff = dobrushin(P, N)
        % DOBRUSHIN:
        % See Section 4.3.2. of Stachurski (2008).
        % If P = obj.P is a Markov matrix, DCOEFF is Dobrushin coefficient.
        % If DCOEFF ~= 0, then P has a (unique) ergodic distribution.
        % If DCOEFF == 0, P may not have an ergodic distribution. But it is
        % possible that DOBRUSHIN(P^N) > 0 for some positive integer N.
            p = P^N;
            ns = size(p,1);
            ergcoeffs = zeros(ns,1);
            for i=1:ns-1
                for k=i+1:size(p,1)
                    ergcoeffs(i) = sum(abs(p(i,:)-p(k,:)));
                end
            end
            dcoeff = 1 - 0.5*max(ergcoeffs);
        end

        function [pi_inf, dcoeff] = ergodist(self,p,t,Ndob)
        % ERGODIST:
        % [PI_INF, DCOEFF] = ERGODIST(P) returns the ergodic distribution
        % and the Dobrushin coefficient.
        %
        % INPUT:
        %       obj     :    class containing P, Markov matrix
        %       p       :    Markov matrix
        %       t       :    optional positive integer, check P is Markov
        %                    matrix
        %       Ndob    :    Some large positive integer, Dobrushin test
        %
        % OUTPUT:
        %       pi_inf  :    ergodic distribution, (1 x size(P,1))
        %       dcoeff  :    Dobrushin coefficient

            ns = size(p,1);

            if ns ~= self.NS
               warning('ERGODIST:Error','size(P) must equal N+M')
            end

            if ~exist('t','var')
                test = 1;
            else
                test = t;
            end
            if test
                if sum(sum(p<0)) ~= 0
                    warning('ERGODIST:Error',...
                                           'Some element in P is negative')
                end
                if sum(sum(p,2)) ~= size(p,1)
                    warning('ERGODIST:Error',...
                                           'Rows in P must be prob. dist.')
                end
                if ~exist('Ndob','var')
                    n = 1;
                else
                    n = Ndob;
                end
                dcoeff = dobrushin(p, n);
                if dcoeff == 0
                    warning('ERGODIST:Error',...
                                              'Dobrushin coefficient zero')
                end
            end

            z = zeros(1,ns);
            z(end) = 1;
            PMI = p - eye(ns);
            PMI(:,end) = ones(ns,1);

            pi_inf = z / PMI;
        end

        function s = initsurplus(self,lambda,bvec)
            % INITSURPLUS:
            % Given lambda calculate surplus required
            s = lambda * bvec' / (1-self.DELTA);
        end

        function vG = swf(self)
            % SWF:
            % Calculate value of government (social welfare) payoff

            vG = self.Results.lambda * self.Results.v';
        end

    end
end
