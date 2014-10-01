function [xI, yI]=findCube2D(xP,yP)
global X Y
    dist=(X(:)-xP).^2+(Y(:)-yP).^2;
    [mdist, mdistIn]=min(dist);
disp(['closest ',num2str(mdistIn)])
    [yI, xI]= ind2sub(size(X),mdistIn); % closest indexes of point
if mdist==0
    return
end
    xsign=sign(xP-X(mdistIn));
    ysign=sign(yP-Y(mdistIn));

     x1c=sort([xI xI+xsign]);
     x2c=sort([yI yI+ysign]);
     [xI, yI]=meshgrid(x1c,x2c);
end