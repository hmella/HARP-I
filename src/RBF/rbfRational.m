
%     Matlab Radial Basis Function Toolkit (MRBFT)
% 
%     Project homepage:    http://www.scottsarra.org/rbf/rbf.html
%     Contact e-mail:      sarra@marshall.edu
% 
%     Copyright (c) 2016-17 Scott A. Sarra
% 
%     Licensing: MRBFT is under the GNU General Public License ("GPL").
%     
%     GNU General Public License ("GPL") copyright permissions statement:
%  ------------------------------------------------------------------------
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
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% ---------------- function summary ---------------------------------------
% -------------------------------------------------------------------------
%
% class rabRational - implements a rational RBF method from the manuscript
%                     "A rational radial basis function method for accurately 
%                      resolving discontinuities and steep gradients" by
%                      S. A. Sarra and Y. Bai.  Under review, Applied
%                      Numerica Mathematics, 2017.

% rbfRationalCoreOne - Called to setup single approximations involving a
%                      set of centers and a shape parameter.  Returns the 
%                      expansion coefficients of p and q to be used in
%                      interpolation and differentiation. The following
%                      functions call this method: interpolate1d,
%                      interpolate2d, diff1dOne, and diff2dOne
%
%
% rbfRationalCoreMultipleSetUp - Sets up for multiple approximations
%                      involving the same centers and shape parameter.
%                      Returns the inverse and Cholesky factorization of B.
% 
% rbfRationalCoreMultiple - Uses the inverse and Cholesky factorization of
%                      B from rbfRationalCoreMultipleSetUp to calculate the
%                      expansion coefficients of p and q.
% 
% interpolate1d - Uses the expansion coefficients of p and q from rbfRationalCoreOne
%                 to interpolate a function of one variable. 
%
% interpolate2d - Uses the expansion coefficients of p and q from rbfRationalCoreOne
%                 to interpolate a function of two variables. 
% 
% diff1dOne - Uses the expansion coefficients of p and q from rbfRationalCoreOne
%             to differentiate a function of one variable. 
%
% diff1dMultiple - Passes the previously calculated inverse and Cholesky
%                  factorization of B to rbfRationalCoreMultiple to get new
%                  expansion coefficients of p and q of to be used to
%                  differentiate functions of one variable using the same 
%                  centers and shape parameter.  For discretizing space 
%                  derivatives of time-dependent PDEs.
% 
% diff2dOne - Uses the expansion coefficients of p and q from rbfRationalCoreOne
%             to differentiate a function of two variables. 
%
% diff2dMultiple - same ad diff1dMultiple except it is for a function of 2 variables.           
% -------------------------------------------------------------------------                
 


classdef rbfRational

methods(Static)
    
% rbfRationalCoreOne
%
% input
%    B    system matrix
%   fc    function values at the center locations
%   mu    optional MDI regularization parameter
%
% output
%  pAlpha  standard RBF exapansion coefficients of the vector p 
%  qAlpha  standard RBF exapansion coefficients of the vector q
    
% called by: interpolate1d, interpolate2d, rbfRational.diff1dOne 
%            and rbfRational.diff2dOne
 
  function [pAlpha, qAlpha] = rbfRationalCoreOne(B,fc,mu)
        if ~exist('mu','var'), mu = 0;  end
        
        N = length(B);
        if mu>0   % MDI regularization
            B(1:N+1:end) = B(1:N+1:end) + mu;  % B = B + mu*eye(N);
        end
 
        L = chol(B,'lower');     % Cholesky factorization of B
        Bi = L'\( L\eye(N) );    % B inverse
        
        Df = diag(fc);
        K = 1.0/sum(fc.^2);
        S = diag(1./( K*fc.^2 + 1 ))*( K*Df*Bi*Df + Bi );
        
        opts.issym = 1;    % S is a real symmetric matrix
        opts.isreal = 1;
        [q,ew] = eigs(S,1,0,opts);    % smallest ew and corresponding ev
        p = Df*q;
        
        % coefficients of the standard RBF interpolant of the vectors p and q
        pAlpha = L'\( L\p );   % solve linear system B*pAlpha = p
        qAlpha = L'\( L\q );   % solve linear system B*qAlpah = q 
 end
    
 % ------------------------------------------------------------------------
 
 % rbfRationalCoreMultipleSetUp
 %
 % input
 %     B   system matrix
 %    mu   optional MDI regularization parameter
 % output
 %     L   Cholesky factorization of B
 %    Bi   inverse of B
 % 
 % example usage: inviscidBurgersRRBF.m, adv1dRRBF.m
 
 
  function [L, Bi] = rbfRationalCoreMultipleSetUp(B,mu)
      if ~exist('mu','var'), mu = 0;  end

      N = length(B);
      if mu>0   % MDI regularization
            B(1:N+1:end) = B(1:N+1:end) + mu;  % B = B + mu*eye(N);
      end
 
      L = chol(B,'lower');     % Cholesky factorization of B
      Bi = L'\( L\eye(N) );    % B inverse
  end
 
 
 % ------------------------------------------------------------------------
 

% rbfRationalCoreMultiple
%
% input
%    B  system matrix
%   Bi  inverse of B (computed by rbfRationalCoreMultipleSetUp)
%    L  Cholesky factorization of B (computed by rbfRationalCoreMultipleSetUp)
%   fc  function value at the center locations
%
% output
%  pAlpha   standard RBF expansion coefficients of the vector p
%  qAlpha   standard RBF expansion coefficients of the vector q
%
% called by: diff1dMultiple and diff2dMultiple
 

 function [pAlpha, qAlpha] = rbfRationalCoreMultiple(B,Bi,L,fc)
        Df = diag(fc);
        K = 1.0/sum(fc.^2);
        S = diag(1./( K*fc.^2 + 1 ))*( K*Df*Bi*Df + Bi );
       
        opts.issym = 1;    % S is a real symmetric matrix
        opts.isreal = 1;
        [q,ew] = eigs(S,1,0,opts);    % smallest ew and corresponding ev

        p = Df*q;
        
        % coefficients of the standard RBF interpolant of the 
        % vectors p and q
        pAlpha = L'\( L\p );   % solve linear system B*pAlpha = p
        qAlpha = L'\( L\q );   % solve linear system B*qAlpah = q  
 end
 

 % ------------------------------------------------------------------------
    
% interpolate1d - 1d interpolation
%
% inputs
%   phi   a RBF (a subclass of rbfx such as iqx or gax)
%     s   shape parameter
%    fc   function values at centers xc
%    xc   centers
%     x   evaluation points
%    mu   optional MDI regularization parameter
%
% output
%    fa   approximate function values at evaluation points x
%
% example usage: stepFunctionRRBF.m
    
    function fa = interpolate1d(phi,s,fc,xc,x,mu)
        if ~exist('mu','var'), mu = 0;  end
       
        r = rbfx.distanceMatrix1d(xc);
        re = rbfx.distanceMatrix1d(xc,x); 
          
        B = phi.rbf(r,s);     % system matrix
        H = phi.rbf(re,s);    % evaluation matrix
        
        [pAlpha, qAlpha] = rbfRational.rbfRationalCoreOne(B,fc,mu);
        fa = (H*pAlpha)./(H*qAlpha);   
    end   
    
    
% interpolate2d - 2d interpolation
%
% inputs
%   phi   a RBF (a subclass of rbfx such as iqx or gax)
%     s   shape parameter
%    fc   function values at centers xc
%    xc   centers
%    yc   
%     x   evaluation points
%     y   
%    mu   optional MDI regularization parameter
%
% output
%    fa   approximate function values at evaluation points x
%
% example usage: interp2dRRBF.m
%    
 
    function fa = interpolate2d(phi,s,fc,xc,x,yc,y,mu)
        if ~exist('mu','var'), mu = 0;  end
        
        r = rbfx.distanceMatrix2d(xc,yc);
        re = rbfx.distanceMatrix2d(xc,yc,x,y); 
        
        B = phi.rbf(r,s);     % system matrix
        H = phi.rbf(re,s);    % evaluation matrix
        
        [pAlpha, qAlpha] = rbfRational.rbfRationalCoreOne(B,fc,mu);
        fa = (H*pAlpha)./(H*qAlpha);    
    end   




% diff1d - 1d differentiation: 1st and 2nd derivatives
%
% inputs
%
%   phi   a RBF (a subclass of rbfx such as iqx or gax)
%     s   shape parameter
%    fc   function values at centers xc
% second  logical variable, if true a second derivative is also returned
%    xc   centers
%
% output
%    fa1   approximate 1st derivative values at evaluation points x
%    fa2   (optional)
%
% example usage: inviscidBurgersRRBF.m, adv1dRRBF.m
        

    function [fa1, fa2] = diff1dMultiple(fc,B,Bi,L,H,H2)
        N = length(fc);
        [pAlpha, qAlpha] = rbfRational.rbfRationalCoreMultiple(B,Bi,L,fc);

        Bp = B*pAlpha;
        Bq = B*qAlpha;
        Hp = H*pAlpha;
        Hq = H*qAlpha;
        
        fa1 = ( Bq.*Hp - Bp.*Hq )./Bq.^2;

         if exist('H2','var')
           Hpp = H2*pAlpha;
           Hqq = H2*qAlpha;
           fa2 = ( 2*Bp.*Hq.^2 + Bq.^2.*Hpp - Bq.*( 2*Hp.*Hq + Bp.*Hqq ) )./Bq.^3;
        end
        
    end
    

    
    
% diff2d - 2d differentiation: 1st and 2nd derivatives
%
% inputs
%   phi   a RBF (a subclass of rbfx such as iqx or gax)
%     s   shape parameter
%    fc   function values at centers (xc,yc)
% second  logical variable, if true a second derivative is also returned
%    xc   centers
%    yc
%
% output
%    fx   approximate 1st derivative values at evaluation points x
%    fy
%    fxx  (optional)
%    fyy  (optional)

function [fx, fy, fxx, fyy] = diff2dMultiple(fc,B,Bi,L,Hx,Hy,H2,H2y)
        N = length(fc);
        [pAlpha, qAlpha] = rbfRational.rbfRationalCoreMultiple(B,Bi,L,fc);

        Bp = B*pAlpha;
        Bq = B*qAlpha;
        
        Hp = Hx*pAlpha;
        Hq = Hx*qAlpha;
        fx = ( Bq.*Hp - Bp.*Hq )./Bq.^2;
        
        Hpy = Hy*pAlpha;
        Hqy = Hy*qAlpha;
        fy = ( Bq.*Hpy - Bp.*Hqy )./Bq.^2;
        
        if exist('H2','var')
           Hpp = H2*pAlpha;
           Hqq = H2*qAlpha;
           fxx = ( 2*Bp.*Hq.^2 + Bq.^2.*Hpp - Bq.*( 2*Hp.*Hq + Bp.*Hqq ) )./Bq.^3;
           
           Hppy = H2y*pAlpha;
           Hqqy = H2y*qAlpha;
           
           fyy = ( 2*Bp.*Hqy.^2 + Bq.^2.*Hppy - Bq.*( 2*Hpy.*Hqy + Bp.*Hqqy ) )./Bq.^3;
        end
        
end
    
% -------------------------------------------------------------------------

% diff1dOne
%
% input
%
%    phi   radial basis function, e.g. iqx, gax
%      s   shape parameter
%     fc   function value at center locations
%     xc   centers
% second   false -> only 1st derivative; true -> second also
%     mu   optional MDI regularization parameter
%
% example usage: diffRRBF.m

 function [fa1, fa2] = diff1dOne(phi,s,fc,xc,second,mu)
        if ~exist('mu','var'), mu = 0;  end
        
        r = rbfx.distanceMatrix1d(xc);   
        B = phi.rbf(r,s);     % system matrix
        H = phi.D1(r,s,r);    % 1st derivative evaluation matrix
        
        [pAlpha, qAlpha] = rbfRational.rbfRationalCoreOne(B,fc,mu);

        Bp = B*pAlpha;
        Bq = B*qAlpha;
        Hp = H*pAlpha;
        Hq = H*qAlpha;
        
        fa1 = ( Bq.*Hp - Bp.*Hq )./Bq.^2;
        
        if second
           H2 = phi.D2(r,s,r);    % 2nd derivative evaluation matrix
           Hpp = H2*pAlpha;
           Hqq = H2*qAlpha;
           fa2 = ( 2*Bp.*Hq.^2 + Bq.^2.*Hpp - Bq.*( 2*Hp.*Hq + Bp.*Hqq ) )./Bq.^3;
        end
        
    end
    

    
    
% diff2d - 2d differentiation: 1st and 2nd derivatives
%
% inputs
%   phi   a RBF (a subclass of rbfx such as iqx or gax)
%     s   shape parameter
%    fc   function values at centers (xc,yc)
% second  logical variable, if true a second derivative is also returned
%    xc   centers
%    yc
%
% output
%    fx   approximate 1st derivative values at evaluation points x
%    fy
%    fxx  (optional)
%    fyy  (optional)
        
    function [fx, fy, fxx, fyy] = diff2dOne(phi,s,fc,xc,yc,second,mu)
        if ~exist('mu','var'), mu = 0;  end

        [r, rx, ry] = rbfx.distanceMatrix2d(xc,yc);   
        B = phi.rbf(r,s);       % system matrix
        Hx = phi.D1(r,s,rx);    % 1st derivative evaluation matrix
        Hy = phi.D1(r,s,ry);
        
        [pAlpha, qAlpha] = rbfRational.rbfRationalCoreOne(B,fc,mu);

        Bp = B*pAlpha;
        Bq = B*qAlpha;
        
        Hp = Hx*pAlpha;
        Hq = Hx*qAlpha;
        fx = ( Bq.*Hp - Bp.*Hq )./Bq.^2;
        
        Hpy = Hy*pAlpha;
        Hqy = Hy*qAlpha;
        fy = ( Bq.*Hpy - Bp.*Hqy )./Bq.^2;
        
        if second
           H2 = phi.D2(r,s,rx);    % 2nd derivative evaluation matrix
           Hpp = H2*pAlpha;
           Hqq = H2*qAlpha;
           fxx = ( 2*Bp.*Hq.^2 + Bq.^2.*Hpp - Bq.*( 2*Hp.*Hq + Bp.*Hqq ) )./Bq.^3;
           
           H2y = phi.D2(r,s,ry);    % 2nd derivative evaluation matrix
           Hppy = H2y*pAlpha;
           Hqqy = H2y*qAlpha;
           
           fyy = ( 2*Bp.*Hqy.^2 + Bq.^2.*Hppy - Bq.*( 2*Hpy.*Hqy + Bp.*Hqqy ) )./Bq.^3;
        end
        
    end


end % methods
   
end % classdef