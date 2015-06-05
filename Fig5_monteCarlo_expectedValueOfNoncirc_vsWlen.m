% Monte Carlo simulation to demonstrate the effect of Hamming window length
% on the circularity of noise subbands.
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

N=2048; %number of samples to use for noncirc sample statistic

M = 512; %length of filter

Ndur = M+N;     %duration of test signal, last N samples are used

Ntrials = 5000; %number of Monte Carlo trials per condition
winlens = [3072 2048 1024 512 256 128 64 32];
bw=winlens;
mu = zeros(Ntrials,length(bw));
cvar = zeros(Ntrials,length(bw));
hvar = zeros(Ntrials,length(bw));
rho = zeros(Ntrials,length(bw));

cvar_o = zeros(Ntrials,length(bw));
hvar_o = zeros(Ntrials,length(bw));
rho_o = zeros(Ntrials,length(bw));
tic;
for ibw=1:length(winlens)
    tic;
    disp(['Running Monte Carlo with window length ' num2str(bw(ibw)) ', ' num2str(ibw) ' of ' num2str(length(bw)) '...']); tic;
    h = hamming(winlens(ibw));
    h = h/sqrt(sum(h.^2));  %make filter unit-energy
for ntrial=1:Ntrials
    x = randn(1,Ndur)+1i.*randn(1,Ndur);    %create circular Gaussian
    y = filter(h,1,x);  %filter the circular Gaussian
    y = y(M+1:end); %clip off filter transient at beginning
    mu(ntrial,ibw) = mean(y);
    cvar(ntrial,ibw) = sum(y.^2)/N; %compute sample compl. variance
    hvar(ntrial,ibw) = sum(abs(y).^2)/N;    %compute sample Hermitian variance
    rho(ntrial,ibw) = cvar(ntrial,ibw)/hvar(ntrial,ibw);    %compute sample noncirc. coeff.
    
    cvar_o(ntrial,ibw) = sum(x.^2)/N; %compute sample compl. variance
    hvar_o(ntrial,ibw) = sum(abs(x).^2)/N;    %compute sample Hermitian variance
    rho_o(ntrial,ibw) = cvar_o(ntrial,ibw)/hvar_o(ntrial,ibw);    %compute sample noncirc. coeff.
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
herr = errorbar(winlens,mnrho,lowererr,uppererr,'b','LineWidth',2);
axis([0 2064 0 1])
set(gca,'FontSize',16);
set(gca,'XTick',fliplr(winlens));

title('$E |\hat{\rho}|$ vs Hamming window length','interpreter','latex');
ylabel('$E |\hat{\rho}|$','interpreter','latex');
xlabel('Hamming window length (number of samples)','interpreter','latex');
set(gca,'XScale','log');


