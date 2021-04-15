function xe = VariableShapeRBFInterp(x,y,p,s,x0,phi)
    
  % RBFinterp performs a RBF interpolation of the data (x,y)
  % at x0, using phi as RBF and p as smoothing factor  



  % Distance matrices
  r  = rbfx.distanceMatrix2d(x(1,:),x(2,:));
  re = rbfx.distanceMatrix2d(x(1,:),x(2,:),...
             x0(1,:),x0(2,:));

  % Variable shape parameters
  sMin = 0.5*s;
  sMax = 1.5*s;
  opt = 3;

  [sn, sm] = phi.variableShape(sMin,sMax,opt,numel(x(1,:)),numel(x0(1,:)));   

  % RBF evaluation and weights estimation
  B = phi.rbf(r,sn);
  w = rbfx.solve(B,y',p,true);

  % RBF evaluation at x0
  Psi = phi.rbf(re,sm);

  % Interpolated values
  xe = Psi*w;

end