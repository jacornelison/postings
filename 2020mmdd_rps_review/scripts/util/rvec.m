function vec=rvec(x)
% vec=rvec(x)
%
% Function to convert any d input to row vector output
%
% A general problem arrises because man Matlab
% functions operate over only 1 dimension, but often
% one wants them to operate over all dimensions
% (eg. sum, mean, std etc.).
% Can get around this by wrapping the operand in a
% a function that makes it 1D. eg:
% 
% std(rvec(x))
%
% Another case where this is useful is where one
% selects ranges of a multi dimensional array and
% wants the result to come out as a row or column vector
% in a single step. eg:
% 
% rvec(img(:,:,3))

vec=x(:)';
