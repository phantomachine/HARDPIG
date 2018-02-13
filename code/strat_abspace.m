function [ ProfileSetA, ProfileSetB, NA, NB ] = strat_abspace(model,...
                                                   Nagrid, Nbgrid, N_Z,M,N)

%% STRAT_ABSET.M
%
% Make a larger cardinality for action sets A and B.
%
% (c) 2012, 2013, Timothy Kam. Email: tcy.kam__at__gmail.com 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Use subject to GNU LGPL licensing terms. Cite this header and 
% author completely in subsequent re-use and modifications.
% =========================================================================
% $Revision: 5.0.0 $  $Date: 2013/09/11 12:45:20 $


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
    


