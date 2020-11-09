function vec = vecp(M)
%get unique elements of symmetric matrix in a vector  
%Order to match Hamimeche-lewis code expectation, see toy model B posting
  dim = size(M, 1);
  vec = zeros(dim*(dim+1)/2, 1);

  counter = 1;

  for iDiag = 0:dim-1
    vec(counter:counter+dim-iDiag-1) = diag(M, iDiag);
    counter = counter + dim - iDiag;
  end
return

