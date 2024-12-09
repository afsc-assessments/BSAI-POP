

postscript(file="q_profile.ps", horizontal=F, family="Times")

par(oma=c(5,5,5,5), mar=c(5,5,5,5), mfrow=c(2,1), cex.axis=1.4, cex.lab=1.4, las=1)


#par(las=1)

#	plot the survey biomass
plot(q_profile_results$q,q_profile_results$AI_srv_biom, type="l", lwd=2, xlab="AI survey q", 
     ylab = "Relative Neg. log likelihood", ylim=c(0,15), col="green")
lines(q_profile_results$q,q_profile_results$EBS_srv_biom, lwd=2, col="purple")
lines(q_profile_results$q,q_profile_results$obj_fun, lwd=2, col="black", lty=2)
legend(0.8,14,c("Obj. Function", "AI srv biomass", "EBS srv biomass" ),
       text.col=c("black","green","purple"),cex=1.2, bty="n")

plot(q_profile_results$q,q_profile_results$fac, type="l", lwd=2, xlab="AI survey q", 
     ylab = "", ylim=c(0,4), col="maroon")
lines(q_profile_results$q,q_profile_results$flc, col="red", lwd=2)
lines(q_profile_results$q,q_profile_results$AI_srv_ac, lwd=2, col="blue")
lines(q_profile_results$q,q_profile_results$EBS_srv_ac, lwd=2, col="orange")
#lines(q_profile_results$q,q_profile_results$obj_fun, lwd=2, col="black", lty=2)
legend(0.9,4,c("Fishery ac", "Fishery lc", "AI srv ac", "EBS srv ac"),
       text.col=c("maroon","red","blue","orange"),cex=1.2, bty="n")

mtext("Relative Neg. log likelihood", side=2, las=3, line=4.0, cex=1.4)

dev.off()




