% Plot to show the expected value of the circularity coefficient of 
% demodulated white noise.
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

x = 0:.0001:.05;
N = 32:4096;
dp = zeros(length(N),length(x));
for ii = 1:length(N)
    x2p = x.*2.*pi;
    dp(ii,:) = sin(N(ii).*x2p./2)./sin(x2p/2).*exp(1i*x2p.*(N(ii)-1)/2)./N(ii);
end

load cmap2
figure;colormap(cmap2);
imagesc(x.*2*pi,N,abs(dp),[0 1]);

colorbar;
set(gca,'FontSize',16);
axis xy;
xlabel('Demodulation Frequency, $\omega$','interpreter','latex');
ylabel('Integration Time, N','interpreter','latex');
title(' $E |\rho_z|$ ','interpreter','latex');
