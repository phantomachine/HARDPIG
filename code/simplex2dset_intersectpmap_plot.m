function simplex2dset_intersectpmap_plot(xPolyVerts, xPolyRand, ...
                                                        xTriIndex, options)
    
NA = numel(xPolyVerts);    
K = numel(xPolyVerts{1});

% Color map
cmap = colormap(lines);                                                 

% Marker options
marker = {  'x','+','*','v','s','^','o','d',... 
            'o','x','+','*','s','d','v','^','<','>','p','h','.',...
            '+','*','o','x','^','<','h','.','>','p','s','d','v',...
            'o','x','+','*','s','d','v','^','<','>','p','h','.'};

if numel(marker) < K
    error('simplex2dset_intersectpmap_inverse_plot:marker',...
            'You need to increase members of marker cell array')
end
        

if strcmp(options.sort_by_K,'on')
    
    for j = 1:K
        
        figure(j)
        for a = 1:NA       
            
            subplot(NA/2,2,a)
            title(strcat('Intersect P(a)(Q(',int2str(j),')) with D; a =',...
                                                            int2str(a)))
            parts = xPolyVerts{a}{j};
            sims = xPolyRand{a}{j};
            
            I = numel(parts);
            
            legenda = cell(I,1);
            
            % Random vectors in polygons:
            for i = 1:I
                hold on
                
                    k = xTriIndex{a}{j}(i);
                    
                    X = sims{i};
                    
                    plot(X(:,1),X(:,2),marker{i},...
                                            'MarkerEdgeColor',cmap(i+a,:));
                                        
                    
                    
                    
                    legenda{i} = strcat('Q(',int2str(k),')');
                    
                hold off
            end
            legend(legenda,'Location','EastOutside')
            
            % Polygon vertices and edges:
            for i = 1:I
                
                v = parts{i};
                hold on
                plot(v(:,1),v(:,2), 'LineStyle','--',...
                                    'Color',cmap(i+a,:));
            
                hold off
            end
            axis tight
        end
        %hold off
    end
end

if strcmp(options.sort_by_A,'on')
    for a = 1:NA
        
        figure(a)
        for j = 1:K       
            
            subplot(K/4,4,j)
            title(strcat('Intersect P(a)(Q(',int2str(j),')) with D; a =',...
                                                            int2str(a)))
            parts = xPolyVerts{a}{j};
            sims = xPolyRand{a}{j};
            
            I = numel(parts);
            
            legenda = cell(I,1);
            
            % Random vectors in polygons:
            for i = 1:I
                hold on
                
                    k = xTriIndex{a}{j}(i);
                    
                    X = sims{i};
                    
                    plot(X(:,1),X(:,2),marker{i},...
                                            'MarkerEdgeColor',cmap(i+a,:));
                                        
                    
                    
                    
                    legenda{i} = strcat('Q(',int2str(k),')');
                    
                hold off
            end
            legend(legenda, 'Location','South',...
                                'Orientation','horizontal')
            
            % Polygon vertices and edges:
            for i = 1:I
                
                v = parts{i};
                hold on
                plot(v(:,1),v(:,2), 'LineStyle','--',...
                                    'Color',cmap(i+a,:));
            
                hold off
            end
            axis tight
        end
        %hold off
    end
end
    
end