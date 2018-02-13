function [V_MIN, V_MAX, ASTAR] = payoff_bounds(self)

        % upper bound on hyperplane distance: a = 0, cons = WAGE + MBAR
        V_MAX = ucons(self, (self.WAGE + self.MBAR)) - phia(self, 0);  
                                                 % No /(1-DELTA);
                                                 % Normalized current value

        % lower bound on payoffs - solve for upper bound on effort, a_max
            syms a
            a_max = solve(  self.DELTA*exp(-a)*V_MAX ...
                                        - (1 - self.DELTA)*(a^self.PHI) );
            ASTAR = double(a_max);

            clear a a_max

        V_MIN = ucons(self, 0) - phia(self, ASTAR);      % NO /(1-DELTA);
                                                 % Normalized current value
                                               