function simplex2dset_intersectpmap_inverse_plot(dt, xPolyQjRand, ...
                                      xTriIndex,xPolyIndexQjRand, options)
    
    NA = numel(xPolyQjRand);   % Number of a profiles 
    K = numel(xPolyQjRand{1}); % Number of state space partition elements
    
    % Color map
    cmap = colormap(lines);                                                 
    
    % Marker options
    % Marker options
    marker = {  'x','+','*','v','s','^','o','d',... 
            'o','x','+','*','s','d','v','^','<','>','p','h','.',...
            '+','*','o','x','^','<','h','.','>','p','s','d','v',...
            'o','x','+','*','s','d','v','^','<','>','p','h','.'};
        
        if numel(marker) < K
            error('simplex2dset_intersectpmap_inverse_plot:marker',...
                    'You need to increase members of marker cell array')
        end
    
    cycle = [ 1 2 3 1 ];     % Simple walk over triangle vertices

    if strcmp(options.sort_by_K,'on')
    for j = 1:K     
        
        % Current triangle Q(J):
        VertIndex = dt.Triangulation(j,:);
        v = dt.X(VertIndex, :);
        
        v = v(cycle,:);
        
        figure(j)
        for a = 1:NA       
            
            subplot(NA/2,2,a)
            title(strcat('Random vectors in Q(',int2str(j),')','; a =', ...
                                                            int2str(a)))

                X = xPolyQjRand{a}{j};  % Random points in Q(j): (Nsim x 3)

                sims = xPolyIndexQjRand{a}{j};
                I = numel(sims);

                legenda = cell(I,1);
                for i = 1:I
                    hold on

                        % Reverse map: indices in X in Q(j) ...
                        select = xPolyIndexQjRand{a}{j}{i};
                        k = xTriIndex{a}{j}(i);
                        if ~isempty(select) || ~isempty(k)
                           
                           legenda{i} = strcat('Q(',int2str(k),')');                             
                           
                           Xs = X(select,:); % All points in X that went to 
                                           % Q(k(i)) under map P(a) on Q(j)
                           plot(Xs(:,1),Xs(:,2),marker{i},...
                                            'MarkerEdgeColor',cmap(i+a,:));
                           
                           
                                    
                        end
                    hold off
                end
                
                legend(legenda, 'Location','South',...
                                'Orientation','horizontal')
                hold on
                    plot(v(:,1),v(:,2),'ro',v(:,1),v(:,2),'r-');
                hold off
                axis tight           
        end
    end
    end 
    
    if strcmp(options.sort_by_A,'on')
     
        figure(a)
        for j = 1:K     
            
            % Current triangle Q(J):
            VertIndex = dt.Triangulation(j,:);
            v = dt.X(VertIndex, :);
        
            v = v(cycle,:);
            
            subplot(K/4,4,j)
            title(strcat('Random vectors in Q(',int2str(j),')','; a =', ...
                                                            int2str(a)))

                X = xPolyQjRand{a}{j};  % Random points in Q(j): (Nsim x 3)

                sims = xPolyIndexQjRand{a}{j};
                I = numel(sims);

                legenda = cell(I,1);
                for i = 1:I
                    hold on

                        % Reverse map: indices in X in Q(j) ...
                        select = xPolyIndexQjRand{a}{j}{i};
                        k = xTriIndex{a}{j}(i);
                        if ~isempty(select) || ~isempty(k)
                           
                           legenda{i} = strcat('Q(',int2str(k),')');                             
                           
                           Xs = X(select,:); % All points in X that went to 
                                           % Q(k(i)) under map P(a) on Q(j)
                           plot(Xs(:,1),Xs(:,2),marker{i},...
                                            'MarkerEdgeColor',cmap(i+a,:));
                           
                           
                                    
                        end
                    hold off
                end
                
                legend(legenda, 'Location','South',...
                                'Orientation','horizontal')
                hold on
                    plot(v(:,1),v(:,2),'ro',v(:,1),v(:,2),'r-');
                hold off
                axis tight           
        end
    end
    end
end