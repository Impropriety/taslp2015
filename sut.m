function [Ayy, K] = sut(y)
% Computes the strong uncorrelating transform and circularity spectrum of
% a two-column vector y
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
y = y.';

ymu = mean(y,2);
bsxfun(@minus,y,ymu);

N = length(y);

Syy = y*y'./N;      %sample Hermitian covariance of mixed sources
Scyy = y*y.'./N;    %sample complementary covariance of mixed sources

if (det(Syy)<eps)
    % Hermitian matrix is singular; must only be one component
    Ayy=[];
    K=[];
    return
end

[U,K,V]=svd(Syy);
Syy_sqrt = U*diag(sqrt(diag(K)))*(V');

C = Syy_sqrt \ Scyy / (Syy_sqrt.');

[U K V] = svd(C);   %take SVD of coherence matrix
D = (V')*conj(U);   %diagonal matrix for computation of Takagi factorization
F = U*diag(sqrt(diag(D)));

Ayy = (F')/Syy_sqrt;
