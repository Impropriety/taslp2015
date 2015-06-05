% Code to compute the time-frequency circularity spectrum on a noisy 
% utterance
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
clear;close all;
load cmap2.mat;
[xin,fs] = audioread('S_02_07.wav');
x1 = xin';
t = [0:length(x1)-1]./fs;

snr_in = -25;
sensor_distance = 0.05;
speed_of_sound = 344;
angle1 = 90;
angle2 = 30;

delay1 = round(sensor_distance./speed_of_sound*cos(angle1/180*pi)*fs);
delay2 = round(sensor_distance./speed_of_sound*cos(angle2/180*pi)*fs);

wlen = 512; %STFT Hamming window length
alen = 1024; %sample size for SUT computation
astep=1; %SUT processing window hop
overlap = wlen-1; %STFT processing window hop
nfft = 2048; %Number of FFT points
step = wlen-overlap; %STFT window hop

x2=randn(1,length(t));
x2 = x2./std(x2);
x2 = x2.*sqrt(var(x1)./(10^(snr_in/10)));

y1 = [x1 zeros(1,delay1+delay2)] + [x2 zeros(1,delay1+delay2)];
y2 = [zeros(1,delay1) x1 zeros(1,delay2)] + [zeros(1,delay2) x2 zeros(1,delay1)];

y1 = y1(:);
y2 = y2(:);
x1 = x1(:);
x2 = x2(:);
n = 0:(length(x1)-1); n = n(:);

h = hamming(wlen);
h = h./sum(abs(h).^2);
[stft_x1, fc, tx_stft]=swSTFT(x1,fs,wlen,nfft,step,h);
stft_x1 = stft_x1(:,wlen+1:end-wlen);
[stft_y1, fc, ty_stft]=swSTFT(y1,fs,wlen,nfft,step,h);
stft_y1 = stft_y1(:,wlen+1:end-wlen);
[stft_y2, fc, ty_stft]=swSTFT(y2,fs,wlen,nfft,step,h);
stft_y2 = stft_y2(:,wlen+1:end-wlen);

snr_out = 10*log10(mean(x1.^2)/mean(x2.^2));

%preallocate storage for circ coeffs
tinds = [1:astep:size(stft_y1,2)-alen];
tlen = length(tinds);
cy1 = zeros(length(fc),tlen);
cy2 = zeros(length(fc),tlen);
tic;
for ii=1:round(length(fc))-1
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
    
    
    for jj = 1:tlen
        SS = [Ssb1([1:alen]+tinds(jj))' Ssb2([1:alen]+tinds(jj))'];
        [Ayy, K] = sut(SS);
        if ~isempty(K)
            cy1(ii,jj)=K(1,1);
            cy2(ii,jj)=K(2,2);
        end
    end
end
colormap hot;
cmap = colormap;
cmap = flipud(cmap);
close all;
figure;colormap(cmap2);imagesc(ty_stft,fc,cy1,[0 1]);
colorbar;
set(gca,'FontSize',16);
axis xy;
xlabel('Time (s)','Interpreter','latex');
ylabel('Frequency (Hz)','Interpreter','latex');
title('First Recovered Circularity Coefficient Magnitude','Interpreter','latex');

figure;colormap(cmap2);imagesc(ty_stft,fc,cy2,[0 1]);
colorbar;
set(gca,'FontSize',16);
axis xy;
xlabel('Time (s)','Interpreter','latex');
ylabel('Frequency (Hz)','Interpreter','latex');
title('Second Recovered Circularity Coefficient Magnitude','Interpreter','latex');


figure;colormap(cmap);imagesc(tx_stft,fc,10.*log10(abs(stft_x1./max(max(stft_x1))).^2),[-60 0])
colorbar;
set(gca,'FontSize',16);
axis xy;
xlabel('Time (s)','Interpreter','latex');
ylabel('Frequency (Hz)','Interpreter','latex');
title('Spectrogram of speech source','Interpreter','latex');

figure;colormap(cmap);imagesc(ty_stft,fc,10.*log10(abs(stft_y1./max(max(stft_y1))).^2),[-60 0])
colorbar;
set(gca,'FontSize',16);
axis xy;
xlabel('Time (s)','Interpreter','latex');
ylabel('Frequency (Hz)','Interpreter','latex');
title(['Spectrogram of speech and noise mixture, SNR=' int2str(snr_in) ' dB'],'Interpreter','latex');