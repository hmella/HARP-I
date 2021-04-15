function [Xphas, Yphas, Zphas] = RefPhaseSmoothing3D(Xpha,Ypha,Zpha,method,RBFFactor,pspace)

    % Create grid
    Isz = size(Xpha,1:3);
    [X, Y, Z] = meshgrid(1:Isz(2),1:Isz(1),1:Isz(3));

    % Get radial basis functions
    phi = feval(method);
    s = RBFFactor;

    % Smooth phase using RBF interpolations
    tf = ~isnan(Xpha) & ~isnan(Ypha) & ~isnan(Zpha);
    fint = RBFInterp3D([X(tf), Y(tf), Z(tf)]',[Xpha(tf), Ypha(tf), Zpha(tf)]',...
            pspace,s(1),[X(tf), Y(tf), Z(tf)]',phi);
    Xphas = NaN(Isz);
    Yphas = NaN(Isz);
    Zphas = NaN(Isz);
    Xphas(tf) = fint(:,1);
    Yphas(tf) = fint(:,2);
    Zphas(tf) = fint(:,3);

end