%% Generate a texture image
% CMYK colorspace
CMYK = zeros(256,4);
CMYK(:,1) = linspace(7,85,256)';
CMYK(:,2) = linspace(7,20,256)';
% RGB colorspace
RGB = zeros(256,3);
RGB(:,1) = 255*(100-CMYK(:,1)).*(100-CMYK(:,4))./10000;
RGB(:,2) = 255*(100-CMYK(:,2)).*(100-CMYK(:,4))./10000;
RGB(:,3) = 255*(100-CMYK(:,3)).*(100-CMYK(:,4))./10000;
RGB = uint8(RGB);
% Draw image
t = repmat(permute(RGB,[1 3 2]),1,256,1);
imwrite(t,'texture.png');

%% Reversed
imwrite(flipdim(t,1),'texture-ud.png');