% Computation of the circularity spectrum on a synthetic harmonic signal
%
% *************************************************************************
% The MIT License (MIT)
% 
% Copyright (c) 2015 Scott Wisdom and Greg Okopal 
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
% *************************************************************************

clear all;
close all;

N = 2048; %number of samples in noncirc estimate 

fs = 8000;
snr_in = -25;
sensor_distance = 0.05;
speed_of_sound = 344;
angle1 = 90;
angle2 =30;

delay1 = round(sensor_distance./speed_of_sound*cos(angle1/180*pi)*fs);
delay2 = round(sensor_distance./speed_of_sound*cos(angle2/180*pi)*fs);

wlen = 512;
overlap = wlen-1;
nfft = 512;
step = wlen-overlap;
t = 0:1/fs:.5;
x1 = zeros(1,length(t));
x2 = zeros(1,length(t));

ft1 = 250;%target fundamental freq

% harmonic tone:
Nharms = 5;
s=zeros(1,length(t));
for kk=1:Nharms
    s = s + sin(2*pi*kk*ft1.*t);    %tone at f0
end
x1=s;
x1 = x1./std(x1);

x2=randn(1,length(t));
x2 = x2./std(x2);
x2 = x2.*sqrt(var(x1)./(10^(snr_in/10)));

snr_out = 10*log10(mean(x1.^2)/mean(x2.^2));

% delays
y1 = [x1 zeros(1,delay1+delay2)] + [x2 zeros(1,delay1+delay2)];
y2 = [zeros(1,delay1) x1 zeros(1,delay2)] + [zeros(1,delay2) x2 zeros(1,delay1)];

y1 = y1(:);
y2 = y2(:);
x1 = x1(:);
x2 = x2(:);
n = 0:(length(x1)-1); n = n(:);

h = hamming(wlen);
h = h./sum(abs(h).^2);
[stft_y1, fc, ty_stft]=swSTFT(y1,fs,wlen,nfft,step,h);
stft_y1 = stft_y1(:,wlen+1:wlen+N);
[stft_y2, fc, ty_stft]=swSTFT(y2,fs,wlen,nfft,step,h);
stft_y2 = stft_y2(:,wlen+1:wlen+N);
[stft_x1, fc, tx_stft]=swSTFT(x1,fs,wlen,nfft,step,h);
stft_x1 = stft_x1(:,wlen+1:wlen+N);
[stft_x2, fc, tx_stft]=swSTFT(x2,fs,wlen,nfft,step,h);
stft_x2 = stft_x2(:,wlen+1:wlen+N);

%preallocate storage for circ coeffs
cy = zeros(2,length(fc));
cx = zeros(2,length(fc));
cmix = zeros(2,length(fc));

tic;
for ii=1:length(fc)-1
    
    if mod(ii,10)==0
        toc;
        disp(['Processed ' num2str(ii) ' of ' num2str(length(fc)) '...']);
        tic;
    end

    fsb = (ii-1)/nfft;
    Ssb1 = stft_y1(ii,:);
    Ssb2 = stft_y2(ii,:);

    %unwrap bulk phase -- mixed
    Ssb1 = abs(Ssb1).*exp(1i.*(unwrap(angle(Ssb1))-cumsum(2*pi*fsb.*ones(1,length(Ssb1)))));
    Ssb2 = abs(Ssb2).*exp(1i.*(unwrap(angle(Ssb2))-cumsum(2*pi*fsb.*ones(1,length(Ssb2)))));
    
    %unwrap bulk phase -- originals
    Sx1 = stft_x1(ii,:);
    Sx2 = stft_x2(ii,:);
    Ssx1 = abs(Sx1).*exp(1i.*(unwrap(angle(Sx1))-cumsum(2*pi*fsb.*ones(1,length(Sx1)))));
    Ssx2 = abs(Sx2).*exp(1i.*(unwrap(angle(Sx2))-cumsum(2*pi*fsb.*ones(1,length(Sx2)))));
    %calc sample circ coeff of originals
    cx(1,ii) = circ(Ssx1');
    cx(2,ii) = circ(Ssx2');
    cmix(1,ii)=circ(Ssb1');
    cmix(2,ii)= circ(Ssb2');
    
    SS = [Ssb1' Ssb2'];
    [Ayy, K] = sut(SS);
    
    cy(1,ii) = K(1,1);
    cy(2,ii) = K(2,2);     
end
if (length(fc) > 1)
    figure;
    plot(fc,(abs(cx(1,:))),'LineWidth',2);hold on;
    plot(fc,(abs(cx(2,:))),'r--','LineWidth',2);
    set(gca,'FontSize',16);
    title('Source circularity coefficient magnitudes','Interpreter','latex');
    xlabel('Frequency (Hz)','Interpreter','latex');
    ylabel('$|\rho|$','Interpreter','latex','FontSize',20);
    legend({'$|\rho_{y1}|$','$|\rho_{y2}|$'},'Interpreter','latex','FontSize',24,'Location','SouthEast')
    axis([0 2000 0 1])
    figure;
    plot(fc,(abs(cmix(1,:))),'LineWidth',2);hold on;
    plot(fc,(abs(cmix(2,:))),'r--','LineWidth',2);
    set(gca,'FontSize',16);
    title('Mixed circularity coefficient magnitudes','Interpreter','latex');
    xlabel('Frequency (Hz)','Interpreter','latex');
    ylabel('$|\rho|$','Interpreter','latex','FontSize',20);
    legend({'$|\rho_{y1}|$','$|\rho_{y2}|$'},'Interpreter','latex','FontSize',24,'Location','SouthEast')    
    axis([0 2000 0 1])
    figure;
    plot(fc,(abs(cy(1,:))),'LineWidth',2);hold on;
    plot(fc,(abs(cy(2,:))),'r--','LineWidth',2);
    set(gca,'FontSize',16);
    title('Separated circularity coefficient magnitudes','Interpreter','latex');
    xlabel('Frequency (Hz)','Interpreter','latex');
    ylabel('$|\rho|$','Interpreter','latex','FontSize',20);
    legend({'$|\rho_{\hat{x}1}|$','$|\rho_{\hat{x}2}|$'},'Interpreter','latex','FontSize',24,'Location','SouthEast')
    axis([0 2000 0 1])
    
    
    fy = abs(fft(y1,nfft));
    fx1 = abs(fft(x1,nfft));
    fx2 = abs(fft(x2,nfft));
    figure;
    plot(fc,20.*log10(fy(1:length(fc))./max(fy)),'k','LineWidth',2);
    hold on;
    plot(fc,20.*log10(fx2(1:length(fc))./max(fy)),'r--','LineWidth',2);
    plot(fc,20.*log10(fx1(1:length(fc))./max(fy)),'b-.','LineWidth',2);
    set(gca,'FontSize',16);
    title('Spectral magnitude of mixture and components','Interpreter','latex');
    xlabel('Frequency','Interpreter','latex');
    ylabel('Magnitude (dB)','Interpreter','latex');
    legend('Y_1(\omega)','V(\omega)','S(\omega)','Location','SouthEast');
    axis([0 2000 -40 0])
    
inds = 1:length(fc);
inds = inds(inds~=17);
inds = inds(inds~=33);
inds = inds(inds~=49);
inds = inds(inds~=65);
inds = inds(inds~=81);
mean(cy(1,inds))
    
end
