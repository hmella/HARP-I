function [Phi] = RBFInterp(X,Y,gamma)

% size of matrix
s = size(F);

% RBF function
rbf = @(r,gamma) = exp(-r.^2/gamma.^2); 

% build matrix
Phi = zeros(s);
for i=1:s(1)*s(2)
    for j=1:s(1)*s(2)
        % distance
        R = sqrt((X(i)-X(j)).^2 + (Y(i)-Y(j)).^2);
        Phi(i,j) = rbf(R, gamma);
    end
end

return;

end
