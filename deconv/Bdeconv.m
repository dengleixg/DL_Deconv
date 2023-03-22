function [I1, I2] = Bdeconv(I,psf,patchI,mode)

[I1,p1] = deconvblind(I,psf,30);

tmp = I1-I;

hlog = fspecial('log',3,1.5);

if mode == 1
tmp = imfilter(tmp,hlog,'symmetric','conv');
end

if mode == 2
tmp = medfilt2(tmp);
end

I2 = tmp + I;

if mode == 3
patchVar = std(patchI)^2;
DoS = 2.5*patchVar;
I2 = imbilatfilt(I1,DoS,20);
end

end