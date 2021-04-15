function [xe, re] = RBFInterp3D(x,y,p,s,x0,phi)
    
  % RBFinterp performs a RBF interpolation of the data (x,y)
  % at x0, using phi as RBF and p as smoothing factor  

  % Distance matrices
  r  = rbfx.distanceMatrix3d(x(1,:),x(2,:),x(3,:));
  re = rbfx.distanceMatrix3d(x(1,:),x(2,:),x(3,:),...
             x0(1,:),x0(2,:),x0(3,:));

  % RBF evaluation and weights estimation
  B = phi.rbf(r,s);
  w = rbfx.solve(B,y',p,true);

  % RBF evaluation at x0
  Psi = phi.rbf(re,s);

  % Interpolated values
  xe = Psi*w;

end