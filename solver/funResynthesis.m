% Resynthesis the denoised image xstar
% Reconstruct the blurred lines image from the solution $\hat w_r$ of the 
% minimization problem by applying the operator $\mathbf{A}$ and by taking 
% the horizontal inverse Fourier transform
% 
%% If you want to reconstruct with other kernels $g$ and $h$ use this code:

% spreadthin=1;
% hthin=exp(-((-S:S)'/spreadthin).^2/2);
% hthin=hthin/sum(hthin);
% gthin=hthin;
% gthinextend = [zeros(M-S,1) ; gthin ; zeros(M-S,1)];
% gthinfour=fftshift(fft(ifftshift(gthinextend)));
% Gthinfour=diag(gthinfour(M+1:end));
% yremp=conv2(wr*Gthinfour,hthin,'valid');

%% By default we reconstruct with the initial kernels $g$ and $h$:

yremp=opA(wr); % Apply the blur operator on the solution wr
% Extend symetrically and compute the horizontal inverse Fourier transform

yrempsym=[fliplr(conj(yremp(:,2:end))) yremp]; % of size WxH
xemp=zeros(H,W);
% inverse Fourier transform on each line
for n2=1:H
    xemp(n2,:)=fftshift(ifft(ifftshift(yrempsym(n2,:))));
end
% Display the deblurred image obtained

f=figure;
subplot(1,2,1);
imagesc(xstar); colormap gray; axis image;
title('Theoritical');
xlabel('$$x^{\sharp}$$','Interpreter','latex');
set(gca,'xtick',[],'ytick',[]);
subplot(1,2,2);
imagesc(xemp); colormap gray; axis image;
title('Empirical');
xlabel('$$\mathcal{F}_1^{-1} \mathbf{A}w$$','Interpreter','latex');
set(gca,'xtick',[],'ytick',[]);
truesize(f,[300 300]);