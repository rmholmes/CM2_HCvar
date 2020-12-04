function [ylin,b] = cfit(x,y)
xvec = [ones(length(x),1) x x.^2 x.^3];
b = xvec\y;
ylin = xvec*b;
end
