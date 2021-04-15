function [r, rx, ry] = distanceMatrix2dNew(xc,yc,x,y)
  xc = xc(:);   yc = yc(:);
  if nargin==2
      rx = xc - xc';
      ry = yc - yc'; 
       r = sqrt( rx.^2 + ry.^2 ); 
  else
      x = x(:);   y = y(:);
      rx = x - xc';
      ry = y - yc'; 
       r = sqrt( rx.^2 + ry.^2 ); 
   end
end