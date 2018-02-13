.. tute2:

%%%%%%%%%%%%%%%%%%%%
Tutorial 2
%%%%%%%%%%%%%%%%%%%%

Here's another longer **MATLAB** example.


Ted's MATLAB solution
-------------------------

.. code-block:: matlab

    function [dt, DT] = simplex2dset_tripart(K, plot_option)

        % SIMPLEX2DSET_TRIPART.M
        %
        % INPUT:
        %
        % K             :  Number of partition elements of 2D unit simplex.
        %                  K must be 4, 16, 64, ... default is 16.
        %
        % plot_option   :  Options: { 'live', 'final', 'none' }. String variable.
        %
        %
        % OUTPUT:
        %
        % dt            :  Structure array:
        %
        %                  (1) dt.X :   Set of extreme points (posibly common) of
        %                               partition elements, in 2D space.
        %
        %                  (2) dt.Triangulation : Indices $t \in \{1,...,K\}$ to 
        %                                         the set dt.X. E.g. dt.X(t,:)
        %                                         yields partition element #t.
        %                                
        %
        % DT            :  Array equivalent to [ dt.X, 1 - sum(dt.X, 2) ]. Each row
        %                  is a point in 2D unit simplex embedded in 3D.
        %
        % (c) 2012, 2013, Timothy Kam. Email: tcy.kam__at__gmail.com 
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % WARNING: Use subject to GNU LGPL licensing terms. Cite this header and 
        % author completely in subsequent re-use and modifications.
        % =======================================================================
        % $Revision: 5.0.0 $  $Date: 2013/09/11 12:45:20 $ 

        % -------------------------------------------------------------------------
        % 1. PARAMETERS & SETTINGS
        % -------------------------------------------------------------------------

        % Check for correct K input ...
        if rem(K,4) ~= 0
            K = 16;
            warning(['K must be powers of 4, Setting default K = ',K])
        end

        N = log(K)/log(4);
        
        % -------------------------------------------------------------------------
        % 2. PERFORM RECURSIVE PARTITIONING:
        % Partitions are equilateral triangles with vertices in Barycentric 
        % coordinates.
        % -------------------------------------------------------------------------
               
        % Pre-allocation structure:
        T_new = struct([]);

        for n = 1:N
            T_new{n} = zeros(3,3,4^n);   
        end
        
        % Initialize with 3D vertices of 2-simplex:
        T_old = eye(3);
        
        % Main loop for recursive partitioning:
        n = 1;
        while n <= N
            
            for m = 1:4^(n-1)
                T_new{n}(:,:,4*m-3:4*m) = simplex_genpart(T_old(:,:,m));
            end
            
            T_old = T_new{n};

            n = n + 1;
        end

        % Desired K(n) = 4^n number of partitions stored in D{n}, n = 1,2,....
        D = struct([]);
        for n = 1:N
            T_temp = permute(T_new{n}, [2,1,3]);
            D{n} = reshape(T_temp,3,3*4^n)';
        end

        % -------------------------------------------------------------------------
        % 3. DELAUNAY TRIANGULATIONS AND OPTIONAL GRAPHICS
        % -------------------------------------------------------------------------

        if strcmp(plot_option, 'final') == 1 || strcmp(plot_option, 'none')
            n = N;
            if strcmp(plot_option, 'final') == 1
                fprintf('\nPlotting only final partition scheme, K = %i ...',K)
            end
        elseif strcmp(plot_option, 'live') == 1
            n = 1;
            fprintf('\nI am plotting all partition schemes ...')
        else
            return
        end

        while n <= N
        
            if strcmp(plot_option, 'none') ~= 1
              fprintf('\nNow plotting partitions with K = %i elements..\n',4^n)
            end
            
            % Create Delaunay complex for each K(n) = 4^n case:
            DelaunayTriangulate = struct([]);
            DT3D = struct([]);

            [ DelaunayTriangulate{n}, DT3D{n}, K_DT, ic ] = ...
                                                        simplex_delaunay(D, n);

            dt = DelaunayTriangulate{n};
            DT = DT3D{n};
            J = 4^n;   

            if strcmp(plot_option, 'none') ~= 1          
                simplex2dset_partdraw(D, T_new, dt, ic, DT, K_DT, K, J, n);
                fprintf('... DONE plotting, K = %i elements\n', 4^n)
            end
            
           n = n + 1;     
        end
    end
