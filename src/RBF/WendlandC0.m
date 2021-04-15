%     Matlab Radial Basis Function Toolkit (MRBFT)
% 
%     Project homepage:    http://www.scottsarra.org/rbf/rbf.html
%     Contact e-mail:      sarra@marshall.edu
% 
%     Copyright (c) 2016 Scott A. Sarra
% 
%     Licensing: MRBFT is under the GNU General Public License ("GPL").
%     
%     GNU General Public License ("GPL") copyright permissions statement:
%     **************************************************************************
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


classdef WendlandC0 < rbfx
    methods
        function obj = WendlandC0()  % constructor
          obj@rbfx();  % call constructor of the superclass 
        end

        function v = rbf(obj,r,s)
            v = (1-r./s).^2;
            v(r>s) = 0.0;
        end % Wendland's C0


        function d = D1(obj,r,s,x), d = []; end
       
        function d = D2(obj, r, s, x), d = []; end
        
        function d = D3(obj, r, s, x), d = []; end
        
        function d = D4(obj, r, s, x), d = []; end
        
        function d = G(obj, r, s, x, y), d = []; end
        
        function d = L(obj, r, s), d = []; end
        
         % x and y not used but required by the abstract function definition in the superclass
        function d = B(obj, r, s, x, y), d = []; end

     
% D12       
% mixed partial derivative
%      D_{xyy}  d = D12( r, s, x, y ) 
% or   D_{yxx}  d = D12( r, s, y, x) 
% depending on the order of the x and y arguments
       
       function d = D12(obj, r, s, x, y), d = []; end
       
       function d = D22(obj, r, s, x, y), d = []; end
      
   end % methods
end  % class
