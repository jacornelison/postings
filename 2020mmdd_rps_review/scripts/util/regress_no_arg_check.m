function b = regress_no_arg_check(y,X)
%REGRESS Multiple linear regression using least squares.
%   B = REGRESS(Y,X) returns the vector B of regression coefficients in the
%   linear model Y = X*B.  X is an n-by-p design matrix, with rows
%   corresponding to observations and columns to predictor variables.  Y is
%   an n-by-1 vector of response observations.
%
%   GPT: 2011 September 13
%   So far in the pipeline, the function b=regress(y,X) is called an
%   exceptionally large number of times, mostly through the subfunction
%   filter_scans>polysub_scans and to a lesser extent through functions
%   elnod and tweak_field_scans.  Since none of these functions require all
%   of the argument checking or optional output of regress, I have simply
%   copied it and removed a lot of the unused stuff.  Normally it wouldn't
%   be such a big deal, but since we call it millions of times for a week
%   of data, it adds up to about 15 percent of the regression time.  We may
%   consider returning to regress at a later date.

% Remove missing values, if any
wasnan = (isnan(y) | any(isnan(X),2));
if any(wasnan)
   y(wasnan) = [];
   X(wasnan,:) = [];
end

% Use the rank-revealing QR to remove dependent columns of X.
[Q,R,perm] = qr(X,0);

% Compute the LS coefficients, filling in zeros in elements corresponding
% to rows of X that were thrown out.
b = zeros(size(X,2),1);
b(perm) = R \ (Q'*y);