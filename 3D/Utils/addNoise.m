function y = addNoise(y,sig,noise)
%% ADDNOISE Add (Gaussian) noise to y with std sig if noise==Gauss.
% y
%    clean data
% sig
%    std for Gaussian noise
% noise
%    String for type of noise (only Gauss).

switch noise
    case 'Gauss'
        y = y + sig*(randn(size(y)) + 1i*randn(size(y)));
end

end