% generate K lines with parameters(tk,ak,nk) where:
% tk is the angles of the line
% ak is the amplitude of the line
% nk is the horizontal offset of the line

%% Lines are blurred with a horizontal filter 
% g and a vertical one h. The image of size WxH degraded
% by some noiselevel is called y.

global h hadj hfourextend Hs M Gfour S;

S=ceil(spread*4)-1;            % half-length of the discrete filter h
Hs=H+2*S;                      % height after adding S pixels beyond border
M=(W-1)/2;                     % half-width of the image
h=exp(-((-S:S)'/spread).^2/2); % the Gaussian filter h
h=h/sum(h);                    % filter normalization
vstar=zeros((H+2*S),W);        % image of the lines horizontally blurred
g=h;                           % same blur g horizontally and h vertically

%% Compute v according to two relations

for k=1:K
    Angle=t_k(k);           % the angle of line with respect to -y axis
    a_k(k)=a_k(k)/max(h);   % after blurring a_k*h is at most equal to 255
    amplitude=a_k(k);       % horizontal offset p_k from the image's center
    offset=p_k(k);          % may varying between -M and M
    for n2=1:(H+2*S)
        t=(0:W-1)-M+tan(Angle)*((H+1)/2+S-n2)-offset;
        vec=zeros(1,W);
        for n=-S:S         % Dirichlet : sin(pi*(t-n))./sin(pi*(t-n)/W);diric
%             vec=vec+g(n+S+1)*arrayfun(@(x) Dirichlet_Kernel(x,W),t-n);
            vec=vec+g(n+S+1)*arrayfun(@(x) diric(x,W),t-n);
        end
        vstar(n2,:)=vstar(n2,:)+amplitude/cos(Angle)*vec/W;
    end
end

%% display the images v',c',y

if plotData
    f=figure;
    subplot(1,3,1);
    imagesc(vstar); colormap gray; axis image;
    title('$$v^{\sharp}$$','Interpreter','latex');
    set(gca,'xtick',[],'ytick',[]);
    subplot(1,3,2);
    imagesc(xstar); colormap gray; axis image;
    title('$$x^{\sharp}$$','Interpreter','latex');
    set(gca,'xtick',[],'ytick',[]);
    subplot(1,3,3);
    imagesc(y); colormap gray; axis image;
    title('$$y$$','Interpreter','latex');
    set(gca,'xtick',[],'ytick',[]);
    warning off;
    truesize(f,[200 200]);
    tightfig;
end

%% compute the FFT of these filters for a later use

gextend = [zeros(M-S,1); g; zeros(M-S,1)];% size 2M+1=W cause we need 
gfour=fftshift(fft(ifftshift(gextend)));     % W Fourier coefficients ^g[m]

Hm=(Hs-1)/2;
hextend=[zeros(Hm-S,1) ; h ; zeros(Hm-S,1)];
hfourextend=fftshift(fft(ifftshift(hextend)));

Gfour=diag(gfour(M+1:end)); % matrix diag containing the Fourier coefs
hadj=conj(flipud(h));       % filter of the adjoint of convolution by h

%% At this stage, we have all we need to perform the minization
% on the data y which is eht degeraded imgae, in order to recover the 
% solution w=Fs and then the parameters of lines. The exact solution
% we want to recover is

% Theoritical formula for wstar
wstar=zeros((H+2*S),W);
for n2=1:(H+2*S)
    for m=-M:M
        wstar(n2,m+M+1)=sum(a_k./cos(t_k).*exp(1i*2*pi*...
            (tan(t_k).*((H+1)/2+S-n2)-p_k).*m/W));
    end
end

%% Theoritical formula for F_1(v)
vfourtheo=zeros(Hs,W);
for m=1:W
    vfourtheo(:,m)=wstar(:,m).*gfour(m);
end

% Compute empirical F_1(v)
vfouremp=zeros((H+2*S),W);
for n2=1:(H+2*S)
    vfouremp(n2,:) = fftshift(fft(ifftshift(vstar(n2,:)))); % horizontal FFT
end

if plotComp
    f=figure;
    subplot(1,2,1);
    imagesc(abs(vfourtheo)); colormap gray; axis image;
    title('Theoritical');
    xlabel('$$\mathcal{F}_1 v^{\sharp}$$','Interpreter','latex');
    set(gca,'xtick',[],'ytick',[]);
    subplot(1,2,2);
    imagesc(abs(vfouremp)); colormap gray; axis image;
    title('Empirical');
    xlabel('$$\mathcal{F}_1 v^{\sharp}$$','Interpreter','latex');
    set(gca,'xtick',[],'ytick',[]);
    truesize(f,[300 300]);
end
