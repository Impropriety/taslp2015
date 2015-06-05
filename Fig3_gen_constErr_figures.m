% Monte Carlo simulation to demonstrate the effect of demodulator mistuning
% on the circularity coefficient
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

set(0,'defaulttextinterpreter','latex');

flag_writeFigs=1;

fs = 8e3;%44.1e3;  %sampling frequency
T = 0.05;      %length in seconds
L = floor(T*fs);    %length in samples
n=0:L-1;

Nsamps = 100;       %number of samples in the phase histogram

ferr = linspace(-30,30,201); %errors in demodulation freq. in Hz
Ntrials = 1000;     %number of Monte Carlo trials per ferr

fc=1e3;
params = { {NaN,fc,0.25,pi/3,L,fs}};

for pp=1:length(params)
    
    fm = params{pp}{1};
    
    Hmn1 = zeros(length(ferr),Ntrials);
    Hstar = zeros(1,Ntrials);
    Rmn1 = zeros(length(ferr),Ntrials);
    cvar = zeros(length(ferr),Ntrials);
    hvar = zeros(length(ferr),Ntrials);
    tic;
    for nn=1:Ntrials
        if mod(nn,floor(Ntrials/10))==0
            toc;
            disp(['Performed ' num2str(nn) ' trials of ' num2str(Ntrials) ' total...']);
            tic;
        end
        [m1 x1a rhoi_star hlp] = gen_improper_modulator_function(params{pp}{:});
        Hstar(nn) = entropy1( hist(angle(m1)+pi,0:(2*pi/Nsamps):2*pi*(Nsamps-1)/Nsamps)./L );
        for ii=1:length(ferr)
            mn1nn = x1a.*exp(-j*2*pi.*n.*(fc+ferr(ii))/fs);
            % calculate sample noncirc. coeff.:
            Rmn1(ii,nn) = sum( mn1nn.^2 )/sum( abs(mn1nn).^2 );
            cvar(ii,nn) = sum( mn1nn.^2 );
            hvar(ii,nn) = sum( abs(mn1nn).^2 );
     end
    end

    pmag = abs(rhoi_star);
    
    D = dirichlet_kernel(L,1.*ferr./(fs/2));    %generate Dirichlet kernel to use for analytic expected value
    
    if isnan(fm)
        Dbs = abs(rhoi_star).*abs(D);
    else
        Dbs_numer = zeros(length(D),length(hlp));
        for mm=0:(length(hlp)-1)
            Dbs_numer(:,mm+1) = ( hlp(mm+1)^2 ).*exp(j*2*pi*mm*2.*ferr).*dirichlet_kernel(L-mm,2.*ferr/fs);
        end
        Dbs_denom = sum(hlp.^2);
        Dbs = abs(rhoi_star).*abs( sum(abs(Dbs_numer),2)./Dbs_denom);
    end

    Rmn1all{pp} = Rmn1; %save off data
    rhoi_star_all{pp} = rhoi_star;
    ferr_all{pp} = ferr;
    fm_all{pp} = fm;
    L_all{pp} = L;
end

figure;
barlocs = [19,53,75,85,101,117,127,149,183];
xax = ferr./fs;
mnRmn1 = abs(mean(Rmn1'));
stdRmn1 = abs(std(Rmn1'));
err = nan(1,length(stdRmn1));
err(barlocs) = stdRmn1(barlocs);
herr = errorbar(xax,mnRmn1,err,'b','LineWidth',2);
hold on;
Dbsdashed = Dbs;
Dbsdashed(5:3:end) = NaN;
ph1=plot(xax,Dbsdashed,'k--');
set(ph1,'LineWidth',3);
ph2=plot([xax(1) xax(end)],[pmag pmag],'r');
set(ph2,'LineWidth',2);
lh=legend([ph1 herr ph2],...
        '$|E\hat{\rho}|$ Analytic   ',...
        '$|\hat{\mu}_{\hat{\rho}}|$ Monte Carlo   ',...
        '$|\rho|$ Actual   ',...
        'Location','NorthEast');
set(lh,'Interpreter','latex');
xlabel('Normalized demodulation error, $\omega_\epsilon$','Interpreter','latex');
    figureHandle = gcf;
    %# make all text in the figure to size 14 and bold
    set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')
    set(gca,'FontSize',14)
lhpos = get(lh,'position');
set(lh,'position',[lhpos(1) lhpos(2)+.05 lhpos(3)+.05 lhpos(4)])
