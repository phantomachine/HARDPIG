function [xcenter, ycenter] = simplex2dset_centroid(T)
    xy_center = sum(T,2)/3;
    xcenter = xy_center(1);
    ycenter = xy_center(2);
end