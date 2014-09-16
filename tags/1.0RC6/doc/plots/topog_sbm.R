# make plot of soil depth function for different values of M
thetaS = 0.6;
Zi = seq(0,1000,5)
Ksat = 100;
MM = c(seq(20,300, by=50),350,500,800,1500);


pdf(file="M_parameter06.pdf")
f = thetaS/MM[1]
Ks = Ksat * exp(-f * Zi)
plot (Ks,-Zi,type="l")
#text(Ks[length(Ks)],-Zi[i],paste(MM[i]))
for (i in 2:length(MM)){
	print(i)
	f = thetaS/MM[i]
	Ks = Ksat * exp(-f * Zi)
	lines (Ks,-Zi)
	text(Ks[length(Ks)* 0.5],-Zi[length(Zi)* 0.5],paste(MM[i]))
}

dev.off()

