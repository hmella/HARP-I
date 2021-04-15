function B = Downsample(A,dims,factor,shift)

    % Default shift is zero
    if nargin < 4
        shift = 0;
    end

    % Number of dimensions of the array
    Ndims = numel(size(A));

    % Array dimensions
    a_dims = 1:Ndims;
    a_dims(dims) = [];

    % Indices
    idx = cell([1 numel(dims)]);
    for i=dims
        idx{i} = (1+shift):factor:size(A,i);
    end
    for i=a_dims
        idx{i} = ':';
    end

    % Donwsample array
    S.subs = idx;
    S.type = '()';
    B = subsref(A,S);

end
