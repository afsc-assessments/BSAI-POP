
####  read into R, make plots, and get model summaries

# get current directory
master="C:/Users/paul.spencer/Work/stock_assess/pop/pop24/Nov model/m_24_2_reweighted/retro"

results_dir <- file.path(master,"results")
setwd(results_dir)

modelnames <- c("m_24_2_r2024","m_24_2_r2023","m_24_2_r2022","m_24_2_r2021","m_24_2_r2020","m_24_2_r2019",
                "m_24_2_r2018","m_24_2_r2017","m_24_2_r2016","m_24_2_r2015","m_24_2_r2014")
                
## read in the *.rdat files
for (i in 1:length(modelnames))
{
         eval(parse(text = paste(modelnames[i]," <- dget('",modelnames[i],".rdat')",sep="")))
}

         
totbiom <- as.data.frame(matrix(data=-9,ncol=length(modelnames),nrow=length(m_24_2_r2022$t.series$year)))
names(totbiom) <- modelnames

spbiom <- as.data.frame(matrix(data=-9,ncol=length(modelnames),nrow=length(m_24_2_r2022$t.series$year)))
names(spbiom) <- modelnames

rec <- as.data.frame(matrix(data=-9,ncol=length(modelnames),nrow=length(m_24_2_r2022$t.series$year)))
names(rec) <- modelnames


for (i in 1:length(modelnames))
{
         tmp <- eval(parse(text = paste(modelnames[i],"$t.series$totbiom")))
         totbiom[1:length(tmp),i] <- tmp
         
          tmp <- eval(parse(text = paste(modelnames[i],"$t.series$spbiom")))
          spbiom[1:length(tmp),i] <- tmp

          tmp <- eval(parse(text = paste(modelnames[i],"$t.series$a3recs")))
          rec[1:length(tmp),i] <- tmp

}

setwd(master)

write.csv(rec,"retro_rec.csv")
write.csv(spbiom,"retro_ssb.csv")


         

