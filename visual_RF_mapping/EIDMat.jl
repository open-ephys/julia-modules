
function EIDMat(rp,isamps,f,sigL,sc,seps)

N       = size(isamps,2);
D       = size(isamps,1);
Ns      = size(rp,2);

blas_set_num_threads(1);
Kxx     = kernelD(isamps,isamps,sigL,sc) + seps*eye(N);
(u,s,v) = svd(Kxx);
s         = diagm(1./max(1e-6,s));
Kxxi    = u*s*v';

#Kxxi    = inv(Kxx);
Kxxf    = Kxxi * f;

Kpx     = kernelD(rp,isamps,sigL,sc) ;
S         = sc - sum(Kpx.*( Kpx * Kxxi), 2);

fx      = Kpx * Kxxf;

S       = sqrt(S);
#g       = (fx - fbest)./S ;

#cdf     = .5 * (1+erf(g/sqrt(2)));
#pdf     = exp(-g.*g/2)/sqrt(2*pi);

#Pf      = S .* (g.*cdf + pdf);

Pf         = fx + 4*S;

return Pf;

end