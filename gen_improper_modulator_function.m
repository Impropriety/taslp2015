function [m1 x1a rhoi_star lpb] = gen_improper_modulator_function(fm,fc,rdi,phi,L,fs)
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
n = (0:(L-1));
N  = 1024;   %order of LPF

if ~isnan(fm)
    lpb = fir1(N,fm/(fs/2));        %filter for signal
    lpb = lpb./sqrt(sum(lpb.^2));
else
    lpb = [];
end

rvar = rdi^2;
ivar = 1;

rhoi_star = (rvar-ivar)/(rvar+ivar)*exp(j*phi)*sign(rvar-ivar); %true noncircularity coeff

mn1 = gen_complex_gaussian(rvar,ivar,phi,L);
if ~isnan(fm)
    m1 = fftfilt(lpb,mn1);
else
    m1 = mn1;
end

fc_off = 0;     %range of support for phase jitter (+/- for uniform dist)
fc_jit_fc = 10;  %cutoff of LPF in Hertz; lower numbers give smoother jitter
lpb_jit = fir1(2*N,fc_jit_fc/(fs/2));        %filter for signal
fc_jito = fc_off.*2.*(rand(1,L)-0.5);
if fc_jit_fc>0
    fc_jit = fftfilt(lpb_jit,fc_jito);
else
    fc_jit = fc_jito;
end

cstar = exp(j*2*pi.*n.*(fc+fc_jit)/fs);
x1a = cstar.*m1;
x1 = real(x1a);

rho1 = (sum(x1a.^2))/sum(abs(x1a).^2);

end