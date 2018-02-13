function H = subgradient(np,L)

% % INPUT:
% % 	* np : Dimensionality of hypersphere S^(np-1)
% %     * L  : Finite elements on S^(np-1)
% % 
% % OUTPUT:
% %     * H : is a finite subset of the hypersphere S^(np-1)
% %
% % DEPENDENCIES:
% %
% %     * EQSP MATLAB Toolbox ( Paul Leopardi, 
% %                             http://sourceforge.net/projects/eqsp/ )
% %
% % Email: mortheus__at__gmail.com
% % =======================================================================
% % $Revision: 4.0.3 $  $Date: 2011/05/04 11:38:20 $
% %
% % See also EQ_POINT_SET

    extra_offset = 1;                           % Option to minimize energy

    H = (eq_point_set(np-1, L, extra_offset))'; % H is finite subset of the 
                                                % np-dimensional sphere: 
                                                % S^(np-1). So, for each i,
                                                % sum( H(i,:).^2, 2 ) = 1.