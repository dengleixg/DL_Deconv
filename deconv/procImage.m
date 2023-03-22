clear
close all
clc

I = dicomread('ST_650.dcm');
figure;imshow(I,[]);title('Original Image');
text(size(I,2),size(I,1)+15, ...
    'Image courtesy of Massachusetts Institute of Technology', ...
    'FontSize',7,'HorizontalAlignment','right');

%crop out the brim to avoid the Gibbs Resonance
Icrop = I(121:574,121:574);


I = double(Icrop);
% I = double(I);
minI = min(I,[],'all'); 
maxI = max(I,[],'all'); 
I = (I-minI)/(maxI-minI);
figure;imshow(I,[]);title('Croped Image');
