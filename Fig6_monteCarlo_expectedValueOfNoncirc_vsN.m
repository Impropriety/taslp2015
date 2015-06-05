% Monte Carlo simulation to demonstrate the effect of sample size on the
% sample circularity coefficent of noise subbands
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

rng(0);

N=1024; %number of samples to use for noncirc sample statistic

M = 512; %length of filter

Ndur = M+N;     %duration of test signal, last N samples are used

Ntrials = 5000; %number of Monte Carlo trials per condition
% bw = [0.005:0.005:0.1]; %list of bandwidths
Nvec = [16 32 64 128 256 512 1024 2048 4096];
mu = zeros(Ntrials,length(Nvec));
cvar = zeros(Ntrials,length(Nvec));
hvar = zeros(Ntrials,length(Nvec));
rho = zeros(Ntrials,length(Nvec));

for iN=1:length(Nvec)
    N = Nvec(iN);
    Ndur = M+N;
    disp(['Running Monte Carlo for N=' num2str(N) ', ' num2str(iN) ' of ' num2str(length(Nvec)) '...']); tic;
    h = hamming(512);
    h = h/sqrt(sum(h.^2));  %make filter unit-energy
for ntrial=1:Ntrials
    x = randn(1,Ndur)+1i.*randn(1,Ndur);    %create circular Gaussian
    y = filter(h,1,x);  %filter the circular Gaussian
    y = y(M+1:M+N); %clip off filter transient at beginning
    mu(ntrial,iN) = mean(y);
    cvar(ntrial,iN) = sum(y.^2)/N; %compute sample compl. variance
    hvar(ntrial,iN) = sum(abs(y).^2)/N;    %compute sample Hermitian variance
    rho(ntrial,iN) = cvar(ntrial,iN)/hvar(ntrial,iN);    %compute sample noncirc. coeff.
    
end
    toc;
end


figure;
barlocs=[3 5 7];
mnrho = mean(abs(rho));
upperstdrho = std(abs(rho))+mnrho;
upperstdrho(upperstdrho>1)=1;
upperstdrho = upperstdrho-mnrho;
lowerstdrho = std(abs(rho));
uppererr = nan(1,length(upperstdrho));
lowererr = nan(1,length(lowerstdrho));
uppererr(barlocs) = upperstdrho(barlocs);
lowererr(barlocs) = lowerstdrho(barlocs);
herr = errorbar(Nvec,mnrho,lowererr,uppererr,'b','LineWidth',2);
axis([0 4896 -.1 1.1])
set(gca,'FontSize',16);
set(gca,'XTick',Nvec);
title('$E |\hat{\rho}|$ vs Sample Size','interpreter','latex');
ylabel('$E |\hat{\rho}|$','interpreter','latex');
xlabel('Sample Size (N)','interpreter','latex');
set(gca,'XScale','log');

