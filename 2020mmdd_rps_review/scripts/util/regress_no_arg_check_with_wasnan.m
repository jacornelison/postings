function b = regress_no_arg_check_with_wasnan(y,X,wasnan)
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
%
%   GPT: 2013 May 6
%   In filter_scans>polysub_scans, we decide not to use any places where X
%   is NaN and also determine whether any(~isnan(y)) before running.  As
%   such, wasnan=isnan(y) already gets computed.  In regress_no_arg_check,
%   recalculating this takes about 20 percent of the regression time.  This
%   function is the same as regress_no_arg_check but uses the precomputed
%   value of wasnan so that it doesn't get recalculated.
%
%%%%%%%%%%% The extra input must one way or another match this: %%%%%%%%%%%
%%%%%%%%%%%%%%%%% wasnan = (isnan(y) | any(isnan(X),2)); %%%%%%%%%%%%%%%%%%

% Remove missing values, if any
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