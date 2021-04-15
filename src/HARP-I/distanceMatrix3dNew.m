function [r, rx, ry, rz] = distanceMatrix3dNew(xc,yc,zc,x,y,z)
    xc = xc(:);   yc = yc(:);  zc = zc(:);
    if nargin==3
        rx = xc - xc';
        ry = yc - yc';
        rz = zc - zc';
        r = sqrt( rx.^2 + ry.^2 + rz.^2 );
    else
        x = x(:);   y = y(:);  z = z(:);
        rx = x - xc';
        ry = y - yc';
        rz = z - zc';
        r = sqrt( rx.^2 + ry.^2 + rz.^2 ); 
     end
end