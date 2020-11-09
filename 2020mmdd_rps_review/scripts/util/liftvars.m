function liftvars(A)
% liftvars(A)
%
% Assigns all fields of the struct A to variables in the caller's scope with
% the same name.
%
% EXAMPLE
%   A.x = 1;
%   liftvars(A)
%   disp(x)
%

for ff=rvec(fieldnames(A))
  assignin('caller', ff{:}, A.(ff{:}));
end
end

