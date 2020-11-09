function model = rps_get_mod_model_alldk(aparam,eta,A,B3,B1,source,rs)
if ~exist('rs','var')
    rs = true;
end

s = size(A);
model = NaN(s(1),13);
param = aparam(1:7);
for i = 1:s(1)
    if length(aparam)>7
        param(7) = aparam(6+i);
    end
   m = rps_get_mod_model_vectors(param,eta,A(i,:),B3(i,:),B1(i,:),source);
   model(i,1:length(m)) = m;
end

if rs
    model = reshape(model,1,[]);
end