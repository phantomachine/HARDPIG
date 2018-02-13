% SIMPLEX_PARTITION.M
%
% Recursive algorithm to generate a sequence of equal-area, equilateral
% triangular partitions of a 2-simplex. This generates NPARTS = 4^N ( N =
% 0,1,2,...) partitions of the simplex.
%
% (c) 2013, T. Kam (tcy.kam@gmail.com)

warning off all
close all
clear 
clc

% -------------------------------------------------------------------------
% 1. PARAMETERS & SETTINGS
% -------------------------------------------------------------------------

    % Set number of partition elements
    K = 16; % Must be 4, 16, 64, 265, ...
    
    % Plotting options:
    plot_all = 'no';
    plot_final = 'yes';

% -------------------------------------------------------------------------
% 2. PERFORM RECURSIVE PARTITIONING:
% Partitions are equilateral triangles with vertices in Barycentric 
% coordinates.
% -------------------------------------------------------------------------
    
N = log(K)/log(4);

        % Check for correct K input ...
        if rem(2*N,2) ~= 0
            error('\n\t ... Set K as one of: 4, 16, 64, 265, ...');   
        end
        
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
% 3. DELAUNAY TRIANGULATIONS AND GRAPHICS
% -------------------------------------------------------------------------

if strcmp(plot_final, 'yes') == 1 && strcmp(plot_all, 'no') == 1
    n = N;
    disp(['I am plotting only the final partition scheme, K = ', ...
                                                   int2str(4^N), ' ...' ])
elseif strcmp(plot_all, 'yes') == 1 && strcmp(plot_final, 'yes') == 1
    n = 1;
    disp('I am plotting all partition schemes ...')
else
    return
end

while n <= N
    
    fprintf('\n\n Now plotting partitions with K = %i elements', 4^n)
    
    % Create Delaunay complex for each K(n) = 4^n case:
    DelaunayTriangulate = struct([]);
    DT3D = struct([]);
    
    [ DelaunayTriangulate{n}, DT3D{n}, K_DT, ic ] = simplex_delaunay(D, n);
    
    dt = DelaunayTriangulate{n};
    DT = DT3D{n};
    J = 4^n;
    simplex_plot; 
        
    fprintf('\n\n\t\t DONE plotting partitions with K = %i elements\n',4^n)
        
   n = n + 1;     
end

