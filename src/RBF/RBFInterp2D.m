function [xe, re, dxe, dye] = RBFInterp2D(x,y,p,s,x0,phi)
    
  % RBFinterp performs a RBF interpolation of the data (x,y)
  % at x0, using phi as RBF and p as smoothing factor  

  % Distance matrices
  r  = rbfx.distanceMatrix2dNew(x(1,:),x(2,:));
  re = rbfx.distanceMatrix2dNew(x(1,:),x(2,:),...
             x0(1,:),x0(2,:));

  % RBF evaluation and weights estimation
  B = phi.rbf(r,s);
  w = rbfx.solve(B,y',p,true);

  % RBF evaluation at x0
  Psi = phi.rbf(re,s);

  % RBF evaluation of derivatives at x0
  Psi_dx = phi.D1(re,s,x(1,:),x0(1,:));
  Psi_dy = phi.D1(re,s,x(2,:),x0(2,:));  

  % Interpolated values
  xe = Psi*w;
  dxe = Psi_dx*w;
  dye = Psi_dy*w;  

end