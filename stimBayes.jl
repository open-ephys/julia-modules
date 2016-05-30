
function stimBayes(nSpikes, iter)
global fprev;
global xprev;

fprev[iter] = sum(nSpikes);

npix 	= [30;30];
D 		= 2;
LB 		= ones(D,1);
UB 		= npix;

x1 = broadcast(+, zeros(npix[1], npix[2]), collect(1:npix[1]));
x2 = broadcast(+, zeros(npix[1], npix[2]), collect(1:npix[2])');
xt          = zeros(D,length(x1));
xt[1,:]     = x1;
xt[2,:]     = x2;

nInit       = 50;
xk          = makeSobol(LB,UB,nInit);

k = iter;

if k < nInit+1
    xprop        = xk[:,k];
else	
	sigL       = 10;  
	seps       = .5;  
	sc         = 2;  
    xprop       = OptProbBayes(xt,xprev[:,max(1,iter-300):iter],fprev[max(1,iter-300):iter],sigL,sc,seps);
end

xprop = round(Int, xprop);

r         = sum(broadcast(-,xt,xprop).^2,1);
sig1     = 1.5;
sig2     = 2*sig1;
Im1     = exp(-r/(2*sig1^2))/sig1;
Im2     = exp(-r/(2*sig2^2))/sig2;
Im       = Im1-Im2;
Im       = Im/maximum(Im)/2;
Im       = 0.5 - Im;
Im       = reshape(Im,npix[1],npix[2]);

xprev[:,iter+1] = xprop;


return Im;
end