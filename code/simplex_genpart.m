function T_new = simplex_genpart(T_old)
    
% SIMPLEX_GENPART.M: Generate Delaunay Triangulations
% (c) T.Kam
%
% INPUT:
%       T_old   : Structure matrix of 3D vertex coordinates
%
% 
% OUTPUT:
%       T_new   : Structure matrix of 3D vertex coordinates of refined
%                 equilateral traingle partition


    if size(T_old,1) ~= 3 && size(T_old,2) ~=3
        error('GENTRIPART.M: input must be a 3 x 3 matrix')
    else
        % Generate new vertices in one-step equal-area split:
        T_temp = zeros(3,3);

        T_temp(1,:) = 0.5*( T_old(1,:) + T_old(2,:) );
        T_temp(2,:) = 0.5*( T_old(2,:) + T_old(3,:) );
        T_temp(3,:) = 0.5*( T_old(1,:) + T_old(3,:) );

        % Store as new set of triangles
        %T_new = struct([]);
         T_new(:,:,1) =  [   T_old(1,:);
                             T_temp(1,:);
                             T_temp(3,:)  ];

         T_new(:,:,2) =  [   T_temp(1,:);
                             T_old(2,:);
                             T_temp(2,:)  ];

         T_new(:,:,3) = [    T_temp(2,:);
                             T_old(3,:);
                             T_temp(3,:)  ];
         
         T_new(:,:,4) = T_temp;
                     
    end   
end