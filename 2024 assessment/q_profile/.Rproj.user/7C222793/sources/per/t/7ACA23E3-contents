
####  read into R, make plots, and get model summaries

# get current directory
master="C:/Users/paul.spencer/Work/stock_assess/pop/pop24/Nov model/m_24_2_reweighted/q_profile"

results_dir <- file.path(master,"results")
setwd(results_dir)

modelnames <- c("pop24_q05",
                "pop24_q06","pop24_q07","pop24_q08","pop24_q09","pop24_q10",
                "pop24_q11","pop24_q12","pop24_q13","pop24_q14","pop24_q15",
                "pop24_q16","pop24_q17","pop24_q18")
                
## read in the *.rdat files
for (i in 1:length(modelnames))
{
         eval(parse(text = paste(modelnames[i]," <- dget('",modelnames[i],".rdat')",sep="")))
}

## define the likelihood components for the profile
likeprof_names <- c("obj_fun","fish.ac","fish.lc","AI_survey.biom","EBS_survey.biom","AI_survey.sac","EBS_survey.sac")
pretty_names <- c("q","obj_fun","fac","flc","AI_srv_biom","EBS_srv_biom","AI_srv_ac","EBS_srv_ac")



q_vec <- seq(0.5,1.8,0.1)
         
q_profile_results <- as.data.frame(matrix(data=-9,nrow=length(modelnames),ncol=length(likeprof_names)+1))
q_profile_results[,1] <- q_vec

for (i in 1:length(modelnames))
{
  for (j in 1:length(likeprof_names)){
     q_profile_results[i,j+1] <- as.numeric(eval(parse(text = paste0(modelnames[i],"$datalikecomp","['",likeprof_names[j],"']"))))
  } 
}
names(q_profile_results) <- pretty_names

# rescale results so min likelihood is at 0 

for (j in 2:(length(likeprof_names)+1)){
  q_profile_results[,j] <- q_profile_results[,j] - min(q_profile_results[,j])
}

setwd(master)


         

