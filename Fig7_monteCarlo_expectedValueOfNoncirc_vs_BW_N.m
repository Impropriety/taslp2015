% Monte Carlo simulation to demonstrate the effect of varying Hamming
% window length and sample size parameters together
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
load cmap2.mat
rng(0);

N=2048;%1024; %number of samples to use for noncirc sample statistic

M = 512; %length of filter

Ndur = M+N;     %duration of test signal, last N samples are used

Ntrials = 1000; %number of Monte Carlo trials per condition
bwvec = [0.005:0.005:0.25]; %list of bandwidths
winlens = [64:64:2304];
bwvec = winlens;
Nvec = [64:64:2304];%[256 512 1024 2048 4096 8192];
mu = zeros(Ntrials,length(Nvec),length(bwvec));
cvar = zeros(Ntrials,length(Nvec),length(bwvec));
hvar = zeros(Ntrials,length(Nvec),length(bwvec));
rho = zeros(Ntrials,length(Nvec),length(bwvec));
mnrho = zeros(length(Nvec),length(bwvec));
for ibw=1:length(winlens)
    winlen = winlens(ibw);
    for iN=1:length(Nvec)
        N = Nvec(iN);
        Ndur = M+N;
        disp(['Running Monte Carlo for N=' num2str(N) ', ' num2str(iN) ' of ' num2str(length(Nvec)) '...']); tic;
        h = hamming(winlen);
        h = h/sqrt(sum(h.^2));  %make filter unit-energy
        for ntrial=1:Ntrials
            x = randn(1,Ndur)+1i.*randn(1,Ndur);    %create circular Gaussian
            y = filter(h,1,x);  %filter the circular Gaussian
            y = y(M+1:M+N); %clip off filter transient at beginning
            mu(ntrial,iN,ibw) = mean(y);
            cvar(ntrial,iN,ibw) = sum(y.^2)/N; %compute sample compl. variance
            hvar(ntrial,iN,ibw) = sum(abs(y).^2)/N;    %compute sample Hermitian variance
            rho(ntrial,iN,ibw) = cvar(ntrial,iN,ibw)/hvar(ntrial,iN,ibw);    %compute sample noncirc. coeff.
        end
        toc;
        mnrho(iN,ibw) = mean(abs(rho(:,iN,ibw)));
    end
end

figure;colormap(cmap2);
imagesc(winlens,Nvec,mnrho,[0 1]);
colorbar;
set(gca,'FontSize',16);
title('$E |\rho|$','interpreter','latex');
xlabel('Hamming window length (number of samples)','interpreter','latex');
ylabel('Sample Size (N)','interpreter','latex');
axis xy;
