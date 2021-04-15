function B = flatten(A,format)

  if nargin < 2
    format = 'col';
  end

  if strcmp(format,'row')
    B = A(:)';
  elseif strcmp(format,'col')
    B = A(:);
  end

end