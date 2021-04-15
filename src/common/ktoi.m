function [Data] = ktoi(Data, Dims)
% ITOK - 
% 

if nargin<2,
  Data = fftshift(ifftn(ifftshift(Data)));
else
  for d=[1:length(Dims)],
    Data = fftshift(ifft(ifftshift(Data),[],Dims(d)));
  end
end
return
