function output = stepcorrchull(H, C_new, K, plot_option)

%% STEPCORRCHULL.M
%
% Construct sample pure strategies supporting initial (lambda0, v0) given
% inner approximating correspondence W^i := [ G(l,k), c(l,k) ]. Uses Bilinear
% Programming approach.
%
% (c) 2012, 2013, Timothy Kam. Email: tcy.kam__at__gmail.com 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Use subject to GNU LGPL licensing terms. Cite this header and 
% author completely in subsequent re-use and modifications.
% =========================================================================
% $Revision: 5.0.0 $  $Date: 2013/09/11 12:45:20 $

    % Vertices of W(k), fit approximate centroid:
    
    N_Z = size(H,2);
    
    VertSet = cell(K,1);
    VertSetChull = VertSet;
    Vert_chull_ind = VertSet;
    Gsearch = VertSet;
    c_normal = VertSet; 
    Wcentroid = zeros(K,N_Z);
    Wcentroid_hull = Wcentroid;
    
    for k = 1:K
        c_k = C_new(:,k);
        %vert_temp = verts_solve(H,c_k);
        
        % Find Vertices of SSE B(W):
        vert_temp = lcon2vert(H,c_k,[],[],1e-8);
            VertSet{k} = vert_temp;
            
        % Q1: Take convhulln of vert_temp? 
        
        ind_temp = convhulln(vert_temp);
        Vert_chull_ind{k} = ind_temp;
        unique_v = unique(ind_temp);
        vert_chull_temp = vert_temp(unique_v,:);
        
        VertSetChull{k} = vert_chull_temp;
        
        % Q2: Then define directions and levels (G,c_in) ? for new set of
        % LP constraints instead of (H,c_new)? 
        
        [ G, cm, Geq, cmeq ] = vert2lcon(vert_chull_temp);
        
        % Normalized/direction vector that are normal to each face
        G_normal = [ G; Geq ];
        lengthG = sqrt(sum(G.^2, 2));
        
        if sum(lengthG) ~= numel(lengthG)
           lengthG = repmat(lengthG,1,N_Z);
           G_normal = G_normal ./ lengthG;
        end
        
        Gsearch{k} = G_normal;
        
        % Levels for each face
        c_normal{k} = [ cm; cmeq ];
        
        % Q3: If we go with Q2, then we can find CENTROID! 
        
        % Roughly, using center of containing hypersphere:
        [Wcentroid(k,:),~] = sphereFit(vert_chull_temp);
        
        % With some finesse, using the vertex points of convex hull to
        % calculate a weighted center of convex polytope:
        
         
         Wcentroid_hull(k,:) = centroid2( vert_temp(unique_v,:) );
    end
    
    %% Pack results into structure OUTPUT
    
    output.VertSet = VertSet;
    output.VertSetChull = VertSetChull;
    output.Vert_chull_ind = Vert_chull_ind;
    output.Gsearch = Gsearch;
    output.c_normal = c_normal; 
    output.Wcentroid = Wcentroid;
    output.Wcentroid_hull = Wcentroid_hull;
    

%%  Show and tell:
if strcmp(plot_option, 'on') == 1;
    figure('name','Original SSE Extreme Points $(H,c)$')
    
    for k = 1:K
       subplot(K/4,K/4,k)
       title(strcat('k = ',int2str(k)))
       plot3(VertSet{k}(:,1), VertSet{k}(:,2), VertSet{k}(:,3), 'o', ...
                    Wcentroid_hull(k,1), Wcentroid_hull(k,2), ...
                                                Wcentroid_hull(k,3), 'or')
       campos([2 13 10])
       xlabel('w_1')
       ylabel('w_2')
       zlabel('w_3')
    end
    
%%  Show and tell:
    figure('name','Vertices of Convex-hull co(W)')
    
    for k = 1:K
       subplot(K/4,K/4,k)
       title(strcat('k = ',int2str(k)))
       plot3(VertSetChull{k}(:,1), ...
                VertSetChull{k}(:,2), ...
                    VertSetChull{k}(:,3), 'o', ...
                    Wcentroid_hull(k,1), Wcentroid_hull(k,2), ...
                                            Wcentroid_hull(k,3), 'or') 
        campos([2 13 10])
        xlabel('w_1')
        ylabel('w_2')
        zlabel('w_3')
    end
    
    %% Interpolated data plot
    figure('name','Convex-hull co(W)')
    for k = 1:K
        x = VertSet{k}(:,1);
        y = VertSet{k}(:,2);
        z = VertSet{k}(:,3);
        ind_temp = Vert_chull_ind{k};
        
          
            subplot(K/4,K/4,k)
            
            title(strcat('k = ',int2str(k)))
            hold on
                plot3(Wcentroid_hull(k,1), Wcentroid_hull(k,2), ...
                                Wcentroid_hull(k,3),... 
                                               'or', 'MarkerFaceColor','r') 
                h = trisurf(ind_temp,x,y,z, 'edgecolor',[0.5 0.5 0.5]);
                alpha(0.25)
                colormap(bone)
                campos([2 13 10])
                %camlight;
                %lighting gouraud;
            hold off
            box off
            set(h,'FaceLighting','gouraud',...
                        'FaceColor','interp',...
                            'AmbientStrength',0.9)
            light('Position',[1 0 0],'Style','infinite');
            set(gca,'color','none')
            axis tight
            xlabel('w_1')
            ylabel('w_2')
            zlabel('w_3')
    end
    fname = strcat('_figures/','stepcorrhull');
            hgsave(fname)
            print('-depsc',fname)
    
end   
 