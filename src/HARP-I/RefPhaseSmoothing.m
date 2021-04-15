function [Xphas, Yphas] = RefPhaseSmoothing(Xpha,Ypha,method,RBFFactor,pspace)

    % Create grid
    Isz = size(Xpha);
    [X, Y] = meshgrid(1:Isz(2),1:Isz(1));

    % Get radial basis functions
    phi = feval(method);
    s = RBFFactor;

    % Smooth phase using RBF interpolations
    tf = ~isnan(Xpha) & ~isnan(Ypha);
    fint = RBFInterp2D([X(tf), Y(tf)]',[Xpha(tf), Ypha(tf)]',...
            pspace,s(1),[X(tf), Y(tf)]',phi);
    Xphas = NaN(Isz(1:2));
    Yphas = NaN(Isz(1:2));
    Xphas(tf) = fint(:,1);
    Yphas(tf) = fint(:,2);

end