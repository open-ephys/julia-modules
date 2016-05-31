function OptProbBayes(xt,ixp,f,sigL,sc,seps)

f    = (f - mean(f))/std(f);

t=1;
Pf  = EIDMat(xt,ixp,f,sigL,sc,seps); 

(pmax,imax) = findmax(Pf);
xprop = xt[:,imax];

return xprop;

end