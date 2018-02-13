function [] = simplex2dset_partdraw(D, T_new, dt, ic, DT, K_DT, K, J, n)

        
    
    colormap = jet(K);

    figure('name', strcat('Ternary Diagram Partition, K = ',' ',...
                                                            int2str(4^n)))

            hold on
                for k = 1:J
                    t = T_new{n}(:,:,k);
                    patch(t(:,1), t(:,2), t(:,3), colormap(k,:));
                    alpha(0.5);
                end

                plot3(D{n}(:,1), D{n}(:,2), D{n}(:,3), 'o', ...
                                                'MarkerFaceColor', 'r')

            hold off
            xlabel('\lambda_{1}')
            ylabel('\lambda_{2}')
            zlabel('\lambda_{3}')
            grid on 
            view(138, 16)


    figure('name', strcat('Delaunay Triangulations (3D view), K = ', ...
                                                     ' ',  int2str(4^n)))
        hold on
            for k = 1:K_DT
                t = dt.Triangulation(k,:);
                patch(DT(t,1), DT(t,2), DT(t,3), colormap(k,:));
                alpha(0.5);
                text(ic(k,1), ic(k,2), 1 - ic(k,1)- ic(k,2), ...
                                 sprintf('T%d', k), ...
                                   'FontWeight', 'bold', ...
                                   'HorizontalAlignment', 'center', ...
                                   'Color','k');
            end

            plot3(D{n}(:,1), D{n}(:,2), D{n}(:,3), 'o', ...
                                                'MarkerFaceColor', 'r')

        hold off
        xlabel('\lambda_{1}')
        ylabel('\lambda_{2}')
        zlabel('\lambda_{3}')
        grid on 
        view(138, 16)

    figure('name', strcat('Delaunay Triangulations (2D plane), K = ',...
                                                         ' ',int2str(4^n)))
        hold on
            for k = 1:K_DT
                t = dt.Triangulation(k,:);
                patch(DT(t,1), DT(t,2), 'w');
                text(ic(k,1), ic(k,2), sprintf('T%d', k), ...
                                   'FontWeight', 'bold', ...
                                   'HorizontalAlignment', 'center', ...
                                   'Color','k');
            end


        hold off
        xlabel('\lambda_{1}')
        ylabel('\lambda_{2}')
        grid on 