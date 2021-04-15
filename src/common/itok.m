function [Data] = itok(Data, Dims)
% ITOK - 
% 
if nargin<2,
  Data = fftshift(fftn(ifftshift(Data)));
else
  for d=[1:length(Dims)],
    Data = fftshift(fft(ifftshift(Data),[],Dims(d)));
  end
end
return

