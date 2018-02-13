function [D, K] = space_distribution(self)

% %-----------------------------------------------------------------------%
% % space_distribution.m
% %-----------------------------------------------------------------------%
% %
% % This function defines the natural product state space of the model. It
% % computes a discretized Cartesian product state space. 
% %
% % =======================================================================
% %     (c) 2011-- T.Kam and R.Stauber.
% %
% % INPUT:
% %
% %     * self    : instance of model class KSMOD2
% %     * .NS     : Number of finite states = card( {-N,...,-1,1,...,M} )
% %     * .PGRID  : Number of elements of set approximating interval [0,1] 
% %     
% % OUTPUT:
% %     * K       : Number of vectors in D
% %     * D       : Finite state space, (K x NS) 
% %
% % DEPENDENCIES:
% %
% %     * CompEcon Toolbox (Miranda-Fackler)
% %         
% % Email: mortheus__at__gmail.com
% % =======================================================================
% % $Revision: 4.0.3 $  $Date: 2011/03/08 11:38:20 $
% %
% % See also GRIDMAKE
  
   
    % partition probability space of state j in [0,1] 
    gridp = linspace(0,1,self.NGRID_P)'; % grids on lambda(j) dimension
      
    D = gridp;
        k = 2;
        while k <= self.NS
            D = gridmake(D,gridp);
        k = k + 1;
        end
           
    % Rule out vectors that are not probability distributions
    D = D( (sum(D,2) == 1), : );
    
    K = size(D,1);
    
end

        
