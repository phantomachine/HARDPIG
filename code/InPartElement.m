function kss = InPartElement(lambda, dt, DtriVerts)

K = size(dt.X,1) + 1;

% Check where lambda belongs to in simplex partition:
    cycle = [ 1 2 3 1 ];
    
    x = lambda(1);
    y = lambda(2);
    
    in = zeros(K,1);
    on = in;
    
    for k = 1:K
        Q_k = DtriVerts(dt.Triangulation(k,:),:);
        Q_k = Q_k(cycle,:);
        xv = Q_k(:,1);
        yv = Q_k(:,2);
        
        [in(k), on(k)] = inpolygon(x,y,xv,yv);
        
    end
    
    % Find index to partition element Q(k) containing lambda_ramsey_ss_max:
    loc1 = find(in == 1);
    loc2 = find(on == 1);
    
    % kss: index of Q(kss) that constains lambda
    if ~isempty(loc1) && isempty(loc2)
        kss = loc1;   
    elseif isempty(loc1) && ~isempty(loc2)
        kss = loc2;
    elseif ~isempty(loc1) && ~isempty(loc2)
        u = rand(1);
        if u < 0.5
            kss = loc1;
        else
            kss = loc2;
        end
    else
        warning('Empty kss!')
    end