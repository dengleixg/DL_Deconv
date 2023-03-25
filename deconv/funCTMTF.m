function [fx_out, MTF_out, V50, V10, theta_out] = funCTMTF(I,pixelL,HVmode)

plotFlag = 0;
index = 1;
limit1 = 50;%mm /after radon
limit2 = 3.12;%40;%mm /after gaussian filter
limit3 = 1.6;%30;%mm /after dff
limit4 = 10;%mm /for fft
Htheta = 0.3;
maxMTF = 0;

theta = -2:0.05:2;

%% Preprocessing
IMask = I;
IMask(:) = normalize(IMask(:),'range');

BW = imbinarize(rot90(IMask));
[BW,~] = imgradient(BW);

%% calcualte the edge angle via Hough Transformation
% theta = -1:0.01:1;
[H,T,~] = hough(BW,'RhoResolution',0.5,'Theta',theta);
peaks = houghpeaks(H,1);

if HVmode == 1%veritcal edge
    theta = T(peaks(2))/180*pi;
end
if HVmode == 2%horizontal edge
    theta = (T(peaks(2))-90)/180*pi;
    IMask = rot90(I);
end

for theta = (theta - Htheta) : 0.01 : (theta + Htheta)
    
    %% get ESF via line integral
    [R,xp] = radon(IMask,theta);
    xp = xp * pixelL;
    R = R(xp>-limit1&xp<limit1);
    xp = xp(xp>-limit1&xp<limit1);
    
    %% ESF with gaussian filter
    w = 2;
    gaussianFilter = exp(-((-(w-1):(w-1))/(w-1)).^2);
    R = conv(R,gaussianFilter,'full');
    R = R(xp>-limit2&xp<limit2);
    xp = xp(xp>-limit2&xp<limit2);
    if plotFlag
        figure(1);subplot(2,2,1);plot(xp,R);title(theta)
    end
    
    %% get LSF by central-difference
    dR = diff(R);
    dx = xp(2:end);
    if plotFlag
        figure(1);subplot(2,2,2);
        axis([-limit4 limit4 min(dR) max(dR)]);
        plot(dx,dR);title(theta);
    end
    y = dR(dx>limit3|dx<-limit3);
    x = dx(dx>limit3|dx<-limit3);
    % remove the baseline
    [fitresult, ~] = createFit(x, y);
    % dR = dR - mean(y);
    try
        dR = dR - (fitresult.a.*dx + fitresult.b);
    catch
        dR = dR - (fitresult.p1.*dx + fitresult.p2);
    end
    if plotFlag
        figure(1);subplot(2,2,2);hold on; plot(dx,dR);hold off
    end
    
    %% add hann window with a frequency domain
    %sampling rate of 0.05cycle/mm
    dR = dR(dx>-limit4&dx<limit4);
    dx = dx(dx>-limit4&dx<limit4);
    w = hann(length(dR));
    dR = dR .* w;
    if plotFlag
        figure(1);subplot(2,2,3);plot(dx,dR);title(theta);
    end
    
    %% get final MTF curve
    MTF = abs(fft(dR));
    MTF = MTF(1:round(end/2));
    MTF = normalize(MTF,'range');
    fx = (0:(length(MTF)-1)) .* ((1/pixelL)/(length(MTF)*2));
    if plotFlag
        h=figure(1);subplot(2,2,4);plot(fx, MTF);title(theta);
%         saveas(h,[num2str(index,'%.4d'),'.png']);
%         pause(0.01)
        index = index + 1;
    end
    
    if MTF(round(end/2))>maxMTF
        theta_out = theta;
        maxMTF = MTF(round(end/2));
        fx_out = fx;
        MTF_out = MTF;
        %saveas(h,[pname,'.png']);
    end
end

V50 = interp1(MTF_out,fx_out,0.5,'pchip');
V10 = interp1(MTF_out,fx_out,0.1,'pchip');
end