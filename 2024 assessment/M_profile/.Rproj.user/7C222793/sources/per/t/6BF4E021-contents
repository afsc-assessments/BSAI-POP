

postscript(file="M_profile.ps", horizontal=T, family="Times")

par(oma=c(5,5,5,5), mfrow=c(1,1), cex.axis=1.4, cex.lab=1.4, las=1)


#par(las=1)

#	plot the survey biomass
plot(M_profile_results$M,M_profile_results$fac, type="l", lwd=2, xlab="Natural Mortality (M)", 
     ylab = "Relative Neg. log likelihood", ylim=c(0,200), col="maroon")
lines(M_profile_results$M,M_profile_results$flc, col="red", lwd=2)
lines(M_profile_results$M,M_profile_results$AI_srv_biom, lwd=2, col="green")
lines(M_profile_results$M,M_profile_results$EBS_srv_biom, lwd=2, col="purple")
lines(M_profile_results$M,M_profile_results$AI_srv_ac, lwd=2, col="blue")
lines(M_profile_results$M,M_profile_results$EBS_srv_ac, lwd=2, col="orange")
lines(M_profile_results$M,M_profile_results$obj_fun, lwd=2, col="black", lty=2)
legend(0.055,150,c("Obj. Function", "Fishery ac", "Fishery lc", "AI srv biomass", "EBS srv biomass", "AI srv ac", "EBS srv ac"),
       text.col=c("black","maroon","red","green","purple","blue","orange"),cex=1.2, bty="n")

dev.off()




