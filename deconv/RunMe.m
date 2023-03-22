clear
close all
clc

%%
param.a = 460;%mm
param.b = 650;%mm
param.c = 196.7385;%mm
param.d = 278;%mm
param.Dd = 9.2;%mGy

%%

I = dicomread('ST_650.dcm') - 1000;
metadata = dicominfo('ST_650.dcm');
figure;imshow(I,[]);title('Original Image');
pixelL = 0.24;%mm
tanLMTF = 5;%mm
tanLCNI = 10;%mm
normLMTF = 4;%mm
normLCNI = 7;%mm
tanLPMTF = round(tanLMTF/pixelL);
tanLPCNI = round(tanLCNI/pixelL);
normPMTF = round(normLMTF/pixelL);
normPCNI = round(normLCNI/pixelL);
%crop out the brim to avoid the Gibbs Resonance
% Icrop = I(121:574,121:574);

% I = double(Icrop);
I = double(I);
minI = min(I,[],'all'); 
maxI = max(I,[],'all'); 
I = (I-minI)/(maxI-minI);
figure;imshow(I,[]);title('Croped Image');

PSF = fspecial('gaussian',11,0.7);

anckerMTF = [432,324];
anckerCNI = [330,308];
% anckerP = [313,206];

X_MTF = [anckerMTF(1)-normPMTF,anckerMTF(1)+normPMTF];
Y_MTF = [anckerMTF(2),anckerMTF(2)+tanLPMTF];
IroiMTF = I(Y_MTF(1):Y_MTF(2),X_MTF(1):X_MTF(2));
figure;imshow(IroiMTF,[]);title('ROI Image for MTF');

X_CNI = [anckerCNI(1),anckerCNI(1)+tanLPCNI];
Y_CNI = [anckerCNI(2)-normPCNI,anckerCNI(2)+normPCNI];
IroiCNI = I(Y_CNI(1):Y_CNI(2),X_CNI(1):X_CNI(2));
figure;imshow(IroiCNI,[]);title('ROI Image for CNI');

CNI_o = funCTCNI(IroiCNI,pixelL,normLCNI,1)

%%
[fx_o, mtf_o, v50o, v10o, ~] = funCTMTF(IroiMTF,pixelL,'original');
figure(100);plot(fx_o,mtf_o,'r.-');hold on;

ai_o = funCTAI(param, v50o, CNI_o)

%% blind deconv

patchI = [];
patchI = [patchI,reshape(I(87:205,287:410),1,[])];
patchI = [patchI,reshape(I(317:375,318:383),1,[])];

[Ib1, Ib2] = Bdeconv(I,PSF,patchI,2);
figure;imshow(Ib1,[])
title('Deblurring with blind PSF')
figure;imshow(Ib2,[])
title('Deblurring with blind PSF(denoise)')

Ib1_mtf = Ib1(Y_MTF(1):Y_MTF(2),X_MTF(1):X_MTF(2));
[fx_b1, mtf_b1, v50b1, v10b1, ~] = funCTMTF(Ib1_mtf,pixelL);
figure(100);plot(fx_b1,mtf_b1,'g.-')

Ib2__mtf = Ib2(Y_MTF(1):Y_MTF(2),X_MTF(1):X_MTF(2));
[fx_b2, mf_b2, v50b2, v10b2, ~] = funCTMTF(Ib2__mtf,pixelL);
figure(100);plot(fx_b2,mf_b2,'b.-')

Ib1 = Ib1*(maxI-minI) + minI;
Ib2 = Ib2*(maxI-minI) + minI;

dicomwrite(int16(Ib1),'Opt1b.dcm',metadata);
dicomwrite(int16(Ib2),'Opt2b.dcm',metadata);

legend({'Original','Opt1b','Opt2b'})

Ib1_cni = Ib1(Y_CNI(1):Y_CNI(2),X_CNI(1):X_CNI(2));
Ib2_cni = Ib2(Y_CNI(1):Y_CNI(2),X_CNI(1):X_CNI(2));

CNI_b1 = funCTCNI(Ib1_cni,pixelL,normLCNI,1)
CNI_b2 = funCTCNI(Ib2_cni,pixelL,normLCNI,1)

ai_b1 = funCTAI(param, v50b1, CNI_b1)
ai_b2 = funCTAI(param, v50b2, CNI_b2)
%%
% figure(101)
% plot(mean(IroiMTF,1),'r.-');hold on;
% plot(mean(Ib1_mtf,1),'g.-');
% plot(mean(Ib2__mtf,1),'b.-');