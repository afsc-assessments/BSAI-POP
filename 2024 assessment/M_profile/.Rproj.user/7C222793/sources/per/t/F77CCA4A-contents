
####  read into R, make plots, and get model summaries

# get current directory
master="C:/Users/paul.spencer/Work/stock_assess/pop/pop24/Nov model/m_24_2_reweighted/M_profile"

results_dir <- file.path(master,"results")
setwd(results_dir)

modelnames <- c("pop24_M01","pop24_M02","pop24_M03","pop24_M04","pop24_M05",
                "pop24_M06","pop24_M07","pop24_M08","pop24_M09","pop24_M10",
                "pop24_M11","pop24_M12","pop24_M13")
                
## read in the *.rdat files
for (i in 1:length(modelnames))
{
         eval(parse(text = paste(modelnames[i]," <- dget('",modelnames[i],".rdat')",sep="")))
}

## define the likelihood components for the profile
likeprof_names <- c("obj_fun","fish.ac","fish.lc","AI_survey.biom","EBS_survey.biom","AI_survey.sac","EBS_survey.sac")
pretty_names <- c("M","obj_fun","fac","flc","AI_srv_biom","EBS_srv_biom","AI_srv_ac","EBS_srv_ac")



M_vec <- seq(0.01,0.13,0.01)
         
M_profile_results <- as.data.frame(matrix(data=-9,nrow=length(modelnames),ncol=length(likeprof_names)+1))
M_profile_results[,1] <- M_vec 


for (i in 1:length(modelnames))
{
  for (j in 1:length(likeprof_names)){
     M_profile_results[i,j+1] <- as.numeric(eval(parse(text = paste0(modelnames[i],"$datalikecomp","['",likeprof_names[j],"']"))))
  } 
}
names(M_profile_results) <- pretty_names

# rescale results so min likelihood is at 0

for (j in 2:(length(likeprof_names)+1)){
  M_profile_results[,j] <- M_profile_results[,j] - min(M_profile_results[,j])
}

setwd(master)


         

