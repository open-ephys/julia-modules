
function makeSobol(LB,UB,nPts)

#srand(123);

d           = size(LB,1);
xt          = rand(d, nPts);   
xt          = round(repmat(LB, 1, nPts) + xt .* repmat(UB-LB, 1, nPts));

return xt
end