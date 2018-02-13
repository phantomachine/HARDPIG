function [] = plot_verts(stepcorr, Z, w, k)

    % stepcorr: class of objects respresenting W
    % Z: extreme points computed by inner ray method
    % w: current w point
    % k: Current index to partition element k = 1,...,K
    Vertices = stepcorr.VertSet{k};

    if size(Z,2) ~= size(Vertices,2) && size(Z,1) == size(Vertices,2)
        Z = Z';
    else
        error('Check the dimensions of Z again')
    end

    Z = [ Z; Z(1,:)];

    % (x,y,z) coordinates of vertices of convex-valued Step Correspondence, W
    x = Vertices(:,1);
    y = Vertices(:,2);
    z = Vertices(:,3);
    % Indices to these(x,y,z) for each domain partition element k
    ind_temp = stepcorr.Vert_chull_ind{k};
    % Set of centroids of convex slice/set W(k)
    Wcentroid_hull = stepcorr.Wcentroid_hull;

        %subplot(K/4,K/4,k)
    figure
        title(strcat('k = ',int2str(k)))
        hold on
            plot3(Wcentroid_hull(k,1), Wcentroid_hull(k,2), ...
                            Wcentroid_hull(k,3),...
                                           'or', 'MarkerFaceColor','r')
            h = trisurf(ind_temp,x,y,z, 'edgecolor',[0.5 0.5 0.5]);
            alpha(0.25)
            colormap(bone)
            campos([2 13 10])
            plot3(x,y,z,'og', 'MarkerFaceColor','none')
            plot3(Z(:,1), Z(:,2), Z(:,3),'s--b', 'MarkerFaceColor','b')
            plot3(w(1), w(2), w(3),'dm', 'MarkerFaceColor','m')
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
