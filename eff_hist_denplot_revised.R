# multiple histogram plot 

#sub_dir  <- c("./Rda_EA/")
#run_date <- c("EA_19990220to20181231_table.comb")

#run_case <- c("w_6_p_80_3D","w_8_p_160_3D",
#              "w_10_p_240_4D","w_12_p_320_5D",

#              "w_6_p_0_TS_2D","w_8_p_0_3D",
#              "w_10_p_0_TS_4D","w_12_p_0_5D",

#              "w_0_p_80_2D", "w_5_p_80_3D",
#              "w_0_p_100_4D","w_5_p_320_5D")

legend.txt <- c("2D: 6m/s,80mm", "3D: 8m/s,160mm", 
                "2D: 6m/s",      "3D: 8m/s",
                "2D: 80mm",      "3D: 160 mm") 

sub_dir  <- c("./Rda_test/")
#
#run_date <- c("EA_19990220to20181231_table.comb")

#

#run_case <- c("w_6_p_80_2D_pos_050_off_000_cut_0", "w_8_p_160_3D_pos_050_off_000_cut_0",
#              "w_6_p_0_2D_pos_050_off_000_cut_0",  "w_8_p_0_3D_pos_050_off_000_cut_0",
#              "w_8_p_80_2D_pos_050_off_000_cut_0", "w_10_p_240_4D_pos_050_off_000_cut_0")
#NOT RUN{
# run_col <- c(rep("red",4),rep("cyan",4),rep("blue",4))
# run_lty <- c(seq(1,4,1),seq(1,4,1),seq(1,4,1))
#}
run_case <- list.files(path=sub_dir, pattern="*.rda")
run_col <- c(rep("red",2),rep("blue",2),rep("cyan",2))
run_lty <- c(seq(1,2,1), seq(1,2,1), seq(1,2,1)) 
 
#0QC1 mean difference  on LAI
#QC2 measurement uncertainty on LAI  
qc1.set <- 1.0
qc2.set <- 0.2
area.set <- 00
####################################################################################
fun_table.qc <- function(wrk.table, qc1.set=0.5, qc2.set=0.2, area.set=1000) {
 
      #subset the datset with QC and QA score
      #Avariable area
      wrk.table <- subset(wrk.table, ((wrk.table$eve.for.pix.aff)> area.set)  )
      #QC1 mean difference  on LAI
      #QC2 measurement uncertainty on LAI  
      #qc1.set <- 0.5
      #qc2.set <- 0.2
      wrk.table <- subset(wrk.table, (abs(wrk.table$qc1.score)<=qc1.set &  abs(wrk.table$qc2.score)>=qc2.set ))
      return(wrk.table)
      } 
###################################################################################


 x_del=0.0
 #load the first case 
 ylim=c(0,30)
 xlim=c(-1.3,1.3)
   
 plot(x=0,  y=0,
      type="l", col="red",xlim=xlim,ylim=ylim,
      xlab="Effect size (ES)",ylab="Count estimate")
 #add zero line 
 abline(v=0, col="lightgray",lty="dashed")
  
  for (icase in 1:length(run_case)) {
  #load the data 
  in_file <- paste(sub_dir,run_case[icase],sep="")
  print(in_file)
  load(in_file)
 #table quality check
 table.comb <- fun_table.qc(wrk.table=table.comb, qc1.set=qc1.set, qc2.set=qc2.set, area.set=area.set) 
 #remove.na 
 print(table.comb$eff.size.for)
 #get density estimate
  dens=density(table.comb$eff.size.for)
  # Plot y-values scaled by number of observations against x values
  lines(x=dens$x-x_del,cex=1.5,
       y=length(table.comb$eff.size.for)*dens$y*dens$bw,
       type="l", col=run_col[icase], lty=run_lty[icase], lwd=2.5)
  }

#add legend 
grid(col="gray",lwd=1.5)
legend("topright", inset=.05, bg="white",
       title="Search Area", 
       legend= legend.txt,
 #      legend= run_case,
       col=run_col, 
       lty=run_lty,cex=1.0)
