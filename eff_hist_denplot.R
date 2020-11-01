# multiple histogram plot 

#sub_dir  <- c("./Rda_4050_Obs6/")
#run_date <- c("WP_20010110to20181231_table.comb")

#run_case <- c("w_3_p_350_TS_123D","w_3_p_350_TS_125D",
#              "w_3_p_350_TS_127D","w_3_p_350_TS_129D",
#              "w_2_p_350_TS_125D","w_3_p_350_TS_125D",
#              "w_5_p_350_TS_125D","w_7_p_350_TS_125D",
#              "w_5_p_40_TS_125D", "w_5_p_80_TS_125D",
#              "w_5_p_100_TS_125D","w_5_p_200_TS_125D",
#              "w_3_p_100_TS_125D" )
#legend.txt <- c("3m/s,p350-3D", "3m/s,p350-5D", "3m/s,p350-7D", "3m/s,p350-9D", 
#               "2m/s,p350-5D", "3m/s,p350-5D", "5m/s,p350-5D", "7m/s,p350-5D",
#                "5m/s,p40-5D", "5m/s,p80-5D","5m/s,p100-5D","5m/s,p200-5D","3m/s,p100-5D"),

sub_dir  <- c("./Rda_pos_off_ace/")
run_date <- c("WP_20010110to20171231_table.comb")

run_case <- c("w_3_p_100_S_125Dpos_060_off_000", "w_3_p_200_S_125Dpos_060_off_000")
# run_col <- c(rep("red",4),rep("cyan",4),rep("blue",4))
# run_lty <- c(seq(1,4,1),seq(1,4,1),seq(1,4,1))
 
 run_col <- c("red","cyan")
 run_lty <- c(1,2)
 
 x_del=0.0
 #load the first case 
 ylim=c(0,50)
 xlim=c(-1.3,1.3)
   
 plot(x=0,  y=0,
      type="l", col="red",xlim=xlim,ylim=ylim,
      xlab="Effect size (ES)",ylab="Count estimate")
 #add zero line 
 abline(v=0, col="lightgray",lty="dashed")
  
  for (icase in 1:2) {
  #load the data 
  in_file <- paste(sub_dir,run_date,run_case[icase],".rda",sep="")
  print(in_file)
  load(in_file)
  table.comb<-subset(table.comb, table.comb$eve.for.pix.aff>5000 & table.comb$date.mon>3 & table.comb$date.mon<=12)
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
 #      legend= legend.txt,
       legend= run_case,
       col=run_col, 
       lty=run_lty,cex=1.0)

