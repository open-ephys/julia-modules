function kernelD(xp0,yp0,sigL,sc)

D  = size(xp0,1);
N  = size(xp0,2); 
M  = size(yp0,2);

#xp = repmat(xp0,1,M);
#yp = reshape(repmat(yp0,N,1),D,N*M);

r = zeros(N,M);
for d = 1:2
	r     = r + (broadcast(-,xp0[d,:]',yp0[d,:])).^2;
end
r  = sqrt(r);
Kn =  (1 + sqrt(3)*r/sigL) .* exp(-sqrt(3)*r/sigL);
 
#Kn  = exp(-r/(2*sigL^2));

Kn = sc * Kn
return Kn

end