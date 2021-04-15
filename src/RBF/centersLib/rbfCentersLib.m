
% The folder rbfX/centersLib contains text files with the locations of 
% centers and evaluation points for various 2d domains.  The function
% rbfCentersLib loads the centers from files into Matlab arrays.
%
% Evaluation points for a domain are numbered as follows:
%    if the centers are for example ch = 4, then the evaluation points for 
%    the same region are numbered as 401, that is, ch = 401


% rbfCentersLib
%
% input
%
%   ch   domain choice
%
% output
%
%  x, y  center locations
%  bpi   indices of the centers located on the boundary of the domain


function [x,y,bpi] =  rbfCentersLib(ch)

   if ch == 1                           % "ameoba" shaped domain
       x = dlmread('x2508.txt','\n');   
       y = dlmread('y2508.txt','\n');
       x = x(:); y = y(:);
       bpi = 1:210;
   elseif ch == 2                       % tilted 8 with a circlular hole
       x = dlmread('x3061.txt','\n');   
       y = dlmread('y3061.txt','\n');
       x = x(:); y = y(:);
       bpi = 1:2;     % NEED TO FIX THIS
   elseif ch == 3
       x = dlmread('x10204.txt','\n');   
       y = dlmread('y10204.txt','\n');
       x = x(:); y = y(:);
       bpi = 8639:9020;
   elseif ch == 4                        % quarter of a circle centers
       XC = dlmread('X618.txt',' ');  
       x = XC(:,1);  y = XC(:,2);
       x = x(:); y = y(:);
       bpi = 1:118;  
   elseif ch == 401                      % quarter of a circle evaluation points
       XC = dlmread('X930.txt',' ');  
       x = XC(:,1);  y = XC(:,2);
       x = x(:); y = y(:);
       bpi = 1:2;                        % no eval points on the boundary     
   end
   
end