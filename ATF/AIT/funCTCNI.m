function CNI = funCTCNI(I,pixelL,normLCNI,Mode)

plotFlag = 0;

[x,y] = size(I);

if Mode == 1%veritcal edge
if x<=9
    error('ROI region is too small.');
end
P = mean(I,1);
S = std(I,0,1);
k = x;
end

if Mode == 2%horizontal edge
if y<9
    error('ROI region is too small.');
end
P = mean(I,2);
S = std(I,0,2);
k=y;
end

y = P(1:round((normLCNI-1)/pixelL))';
x = 1:round((normLCNI-1)/pixelL);

% [fitresult, ~] = createFit(x, y);
% if strcmp(version('-release'), '2021a')
%         Ptmp = P - (fitresult.b.*(1:length(P)) + fitresult.a)';
% else
%         Ptmp = P - (fitresult.p1.*(1:length(P)) + fitresult.p2)';
% end

Ptmp = P;

dP = zeros(size(Ptmp));
d2P = zeros(size(Ptmp));

for m = 5:k-4
    dP(m) = sum(Ptmp(m:m+4))/5-sum(Ptmp(m-4:m-1))/4;
end

for m = 5:k-5
    d2P(m) = dP(m+1)-dP(m);
end

d2Ptmp = d2P(5+1:k-4-1);

[d2Pmax, imax] = max(d2Ptmp);
[d2Pmin, imin] = min(d2Ptmp);

Pmax = P(d2P == d2Pmax);
Pmin = P(d2P == d2Pmin);
Smax = S(d2P == d2Pmax);
Smin = S(d2P == d2Pmin);

CNI = abs(Pmax-Pmin)/sqrt((Smax^2+Smin^2)/2);

if plotFlag
    figure;
    subplot(2,2,1);plot(P,'r');
    hold on;
    plot(imax+5,Pmax,'o');
    plot(imin+5,Pmin,'o');
    title('P');
    subplot(2,2,2);plot(S,'g');
    hold on 
    plot(imax+5,Smax,'o');
    plot(imin+5,Smin,'o');
    title('S');
    subplot(2,2,3);plot(dP,'b');title('dP');
    subplot(2,2,4);plot(d2P,'c');title('d2P');
end

end