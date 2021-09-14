#This R code gives example to implement major results discussed in the manuscipt.
#The codes of hidden Markov modelling are learned from Huang (2018).
#For the access to the codes of spectrum analysis, please contact the authors Komarzynski (2018).

library(depmixS4) #Package for HMM
library(lubridate) #Package for dates and times
library(scales) #Package for plots
library(ggplot2) #Package for plots

source('oneHarmonic_HMMs.R')
source('one_day_summary.R')
source('twoHarmonic_HMMs.R')
source('test_24h_consecutive.R')
#The codes for hidden Markov modelling are learned from Huang et al.(2018).

##Read files
id <- 1101
info<-read.csv(paste0("subject_",id,".csv")) #Read the 5-min average physical activity and chest surface temperature record
diary<-read.csv(paste0("subject_",id,"_diary_sleep.csv")) #Read diary
time<- as.POSIXct(info$time, format="%Y-%m-%d %H:%M:%S", "GMT")  #Get appropriate datetime class
time_date<-day(time)
work_time<- as.POSIXct(diary$Date.time[which(diary$Diary=='work')],format = "%d/%m/%Y %H:%M","GMT")

###################################################################
###Harmonic HMMs (circadian harmonic over the entire study session) 
###################################################################
sf <- 1/12 #Sampling frequency in hour unit, i.e. here sf=5min/60min=1/12
lag<-0:(nrow(info)-1)
sin_baseline<-sin(2*pi*lag*sf/24)
cos_baseline<-cos(2*pi*lag*sf/24)
y_1hmc<-depmix(sqValue~1,transition=~sin_baseline+cos_baseline,
               data=info,nstates=3,family=gaussian(),ntimes=nrow(info)) #Assume that there are 3 hidden states and that the emission density is Gaussian

AIC <- array(NA,10) #The likelihood may have local maxima, thus try 10 starting values and compare AIC
for(a in 1:10){
  set.seed(a)
  HMMfit_1hmcs<-try(fit(y_1hmc,verbose=F),TRUE)
  if(isTRUE(class(HMMfit_1hmcs)=="try-error")) 
  { next } else {AIC[a]=AIC(HMMfit_1hmcs)} 
}
AIC <- as.data.frame(AIC)
set.seed(which(AIC==min(AIC,na.rm=TRUE))) #Select the one with the lowest AIC
HMMfit_1hmcs<-fit(y_1hmc,verbose=F)
HMM_results_1hmc<-oneHarmonic_HMMs(HMM=HMMfit_1hmcs,sin1=sin_baseline,cos1=cos_baseline) #Summarise useful parameters 

##Compute day profile and circadian parameters
circadian_states_prob<-HMM_results_1hmc$circadian_states_prob 
circadian_states_prob$time<-info$time
transition_prob<-HMM_results_1hmc$transition_prob
hour_day_start_range<- seq(12,20,2)
one_day_summary_subject<-one_day_summary(circadian_states_prob=circadian_states_prob[seq(1, nrow(info)),],
                                         transition_prob=transition_prob[seq(1, nrow(info)),],
                                         hour_day_start_range=hour_day_start_range,sf=sf, figure=F) 
HMM_results_1hmc$one_day_summary_subject<-one_day_summary_subject
save(HMM_results_1hmc,file=paste0("HMC_",id,"_new.Rdata"))

##Produce Figure 3A.a, 3A.c and 3A.d
df <- data.frame(time) 
time2 <- c(time[-1],NA)
df$time2 <- time2
df$prob1<- HMM_results_1hmc$prob_ML_states[,1]
df$prob12<- HMM_results_1hmc$prob_ML_states[,1]+HMM_results_1hmc$prob_ML_states[,2]
df$prob123<- HMM_results_1hmc$prob_ML_states[,1]+HMM_results_1hmc$prob_ML_states[,2]+HMM_results_1hmc$prob_ML_states[,3]
df$sqValue <- info$sqValue
df$temp <- info$Temp_processed #Chest temperature
df$sleep <- rep(0,nrow(df)) #This column is to identify sleep, work, and other time periods recorded in the diary
sleep_time<- as.POSIXct(diary$Date.time[which(diary$Diary=='sleep')],format = "%d/%m/%Y %H:%M","GMT") #Sleep periods in the diary
for (j in 1:length(sleep_time)){
  if (j %%2 ==1){
    findsleep <- df[which((date(df$time) == date(sleep_time[j]))& (hour(df$time)==hour(sleep_time[j]))),]
    for (k in 1:nrow(findsleep)){
      findsleep$difference[k] <- abs(as.numeric(findsleep$time[k]-sleep_time[j]))
    }
    start <- findsleep[which(findsleep$difference==min(findsleep$difference)),]$time
    m <- j+1
    findsleep <- df[which((date(df$time) == date(sleep_time[m]))& (hour(df$time)==hour(sleep_time[m]))) ,]
    if (nrow(findsleep)>0){ 
      for (k in 1:nrow(findsleep)){
        findsleep$difference[k] <- abs(as.numeric(findsleep$time[k]-sleep_time[m]))
      }
      end <- findsleep[which(findsleep$difference==min(findsleep$difference)),]$time
    } else{
      end <- df$time[length(df$time)]
    }
    start <- which(df$time==start)
    end <- which(df$time==end)
    df$sleep[start:end] <- 1
  }
}
work_time<- as.POSIXct(diary$Date.time[which(diary$Diary=='work')],format = "%d/%m/%Y %H:%M", "GMT") #Work periods in the diary
for (j in 1:length(work_time)){
  if (j %%2 ==1){
    findwork <- df[which((date(df$time) == date(work_time[j]))& (hour(df$time)==hour(work_time[j]))) ,]
    if (nrow(findwork)>0){
      for (k in 1:nrow(findwork)){
        findwork$difference[k] <- abs(as.numeric(findwork$time[k]-work_time[j]))
      }
      start_time <- findwork[which(findwork$difference==min(findwork$difference)),]$time
      m <- j+1
      findwork <- df[which((date(df$time) == date(work_time[m]))& (hour(df$time)==hour(work_time[m]))) ,]
      if (nrow(findwork)>0){ 
        for (k in 1:nrow(findwork)){
          findwork$difference[k] <- abs(as.numeric(findwork$time[k]-work_time[m]))
        }
        end_time <- findwork[which(findwork$difference==min(findwork$difference)),]$time
      } else{
        end_time <- df$time[length(df$time)]
      }
      start <- which(df$time==start_time)
      end <- which(df$time==end_time)
    } else{
      m <- j+1
      findwork <- df[which((date(df$time) == date(work_time[m]))& (hour(df$time)==hour(work_time[m]))) ,]
      if (nrow(findwork)>0){
        start_time <- df$time[1]
        for (k in 1:nrow(findwork)){
          findwork$difference[k] <- abs(as.numeric(findwork$time[k]-work_time[m]))
        }
        end_time <- findwork[which(findwork$difference==min(findwork$difference)),]$time
        start <- which(df$time==start_time)
        end <- which(df$time==end_time)
      } else {
        start <- 0
        end <- 0
      } 
    }
    if (end > 0){df$sleep[start:end] <- 2}
  }
}
df$sleep<- as.factor(df$sleep)
df$temp[is.na(df$sqValue)] <- NA
df$state <- HMM_results_1hmc$ML_states
df$state[df$state==1] <- HMM_results_1hmc$obs_density_sq$mean[1] #Data for the yellow/orange line in Figure 3A.a 
df$state[df$state==2] <- HMM_results_1hmc$obs_density_sq$mean[2]
df$state[df$state==3] <- HMM_results_1hmc$obs_density_sq$mean[3]

ggplot(df)+ 
  geom_rect(aes(xmin=time,xmax=time2,ymin=0, ymax=25, fill=days), alpha =0.1)+
  geom_point(aes(x=time,y=sqValue), col = "black")+
  geom_point(aes(x=time,y=(temp-35)*2+35-15), col = "brown")+
  geom_line(aes(x=time,y=state), col = "orange")+
  geom_rect(aes(xmin = time,xmax=time2, ymin = -0.7,ymax=-0.4, fill = sleep)) +
  geom_point(aes(x=df$time[1], y=25), colour="black",size=3)+
  geom_point(aes(x=df$time[500], y=25), colour="brown",size=3)+ 
  annotate("rect", xmin = df$time[900], xmax = df$time[940], ymin = 24.8, ymax = 25,fill = "orange",alpha=1)+ 
  annotate("text", x=df$time[10], y=25, hjust=0,label= "5-min average PA (square root)",size=5)+
  annotate("text", x=df$time[510], y=25, hjust=0,label= "5-min average Chestemp",size=5)+
  annotate("text", x=df$time[950], y=25, hjust=0,label= "Most likely states",size=5)+
  scale_fill_manual(name="Diary:",values=c("grey","blue","red"),
                    labels = c("Other", "Sleep", "Work"))+
  theme_bw()+
  theme(text = element_text(size=25),
        axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0,size=20),
        axis.title.y = element_text(color ="black"),
        axis.text.y = element_text(color = "black",size=20),
        axis.title.y.right = element_text(color="brown",size=20),
        axis.text.y.right = element_text(color="brown",size=20),
        plot.title = element_text(face = "bold"),
        legend.position = "right",
        legend.background = element_rect(color = "black", size = 0.3, linetype = "solid"), 
        legend.direction = "vertical",
        legend.title=element_text(size=15),
        legend.text=element_text(size=15))+
  scale_x_datetime(date_breaks = "12 hour", labels = date_format("%d/%m %H:%M"), name="Time")+
  scale_y_continuous(limits=c(-0.7,25),breaks=seq(0,25,5),name="PA (accelerations)",label=scales::number_format(accuracy = 1),
                     sec.axis=sec_axis(~((.+15-35)/2+35), name="Chestemp (¡ãC)"))+
  ggtitle('a') #Figure 3A.a (15 and 35 are arbitrarily selected to plot two different scales)

ggplot(df, aes(x=time)) + 
  geom_ribbon(data=df, 
              aes(ymin=prob1,ymax=prob12), fill="pink", alpha=0.5)+
  geom_ribbon(data=df, 
              aes(ymin=0,ymax=prob1), fill="blue", alpha=0.5)+
  geom_ribbon(data=df, 
              aes(ymin=prob12,ymax=prob123), fill="red", alpha=0.5)+
  theme_bw()+
  geom_rect(aes(xmin = time,xmax=time2, ymin = -0.05,ymax=-0.01, fill = sleep)) +
  scale_x_datetime(date_breaks = "12 hour", labels = date_format("%d/%m %H:%M"), name="Time")+
  scale_y_continuous(breaks=seq(0,1,0.2),limits=c(-0.05,1.1),name="Probability of every state")+
  theme(text = element_text(size=25),
        legend.position="none",
        axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0,size=20),
        plot.title = element_text(face = "bold"))+
  scale_fill_manual(name="Diary:",values=c("grey","blue","red"),
                    labels = c("NA", "Sleep", "Work"))+
  annotate("rect", xmin = df$time[1], xmax = df$time[90], ymin = 1.03, ymax = 1.05,fill = "blue",alpha=0.5)+
  annotate("rect", xmin = df$time[300], xmax = df$time[390], ymin = 1.03, ymax = 1.05,fill = "pink",alpha=0.7)+
  annotate("rect", xmin = df$time[300], xmax = df$time[390], ymin = 1.08, ymax = 1.1,fill = "red",alpha=0.5)+
  annotate("rect", xmin = df$time[700], xmax = df$time[790], ymin = 1.03, ymax = 1.05,fill = "red")+
  annotate("rect", xmin = df$time[700], xmax = df$time[790], ymin = 1.08, ymax = 1.1,fill = "blue")+
  annotate("text", x=df$time[100], y=1.05, hjust=0,label="Inactive state",size=5)+
  annotate("text", x=df$time[400], y=1.05,hjust=0,label= "Moderately active state",size=5)+
  annotate("text", x=df$time[400], y=1.1, hjust=0,label= "Highly active state",size=5)+
  annotate("text", x=df$time[800], y=1.05,hjust=0,label= "Work in diary",size=5)+
  annotate("text", x=df$time[800], y=1.1, hjust=0,label= "Sleep in diary",size=5)+
  ggtitle('c') #Figure 3A.c

circadian_states_prob_forplot <- circadian_states_prob[1:(288*2),]
start_hour<- 20 #Make sure that the daily profile starts around 20:00
start <- min(which(hour(circadian_states_prob_forplot$time)==start_hour))
end <- start+24/sf-1
one_day_prob <- circadian_states_prob_forplot[start:end,]
one_day_prob$time <- as.POSIXct(one_day_prob$time, format="%Y-%m-%d %H:%M:%S", "GMT")
day(one_day_prob$time)[1:(min(which(hour(one_day_prob$time)==0))-1)] <- min(unique(day(one_day_prob$time)))
if (length(unique(day(one_day_prob$time)))==1){
  day(one_day_prob$time)[min(which(hour(one_day_prob$time)==0)):length(one_day_prob$time)] <- unique(day(one_day_prob$time))+1
} else{
  day(one_day_prob$time)[min(which(hour(one_day_prob$time)==0)):length(one_day_prob$time)] <- max(unique(day(one_day_prob$time)))
}

ggplot(one_day_prob, aes(x=time)) +
  geom_ribbon(aes(ymin=state_1,ymax=state_1+state_2), fill="pink", alpha=0.5)+
  geom_ribbon(aes(ymin=0,ymax=state_1), fill="blue", alpha=0.5)+
  geom_ribbon(aes(ymin=state_1+state_2,ymax=state_1+state_2+state_3), fill="red", alpha=0.5)+
  annotate("rect", xmin = one_day_prob$time[1], xmax = one_day_prob$time[20], ymin = 1.03, ymax = 1.05,fill = "blue",alpha=0.5)+
  annotate("rect", xmin = one_day_prob$time[100], xmax = one_day_prob$time[120], ymin = 1.03, ymax = 1.05,fill = "pink",alpha=0.7)+
  annotate("rect", xmin = one_day_prob$time[100], xmax = one_day_prob$time[120], ymin = 1.08, ymax = 1.1,fill = "red",alpha=0.5)+
  annotate("text", x=one_day_prob$time[25], y=1.05, hjust=0,label="IA state",size=5)+
  annotate("text", x=one_day_prob$time[125], y=1.05,hjust=0,label= "MA state",size=5)+
  annotate("text", x=one_day_prob$time[125], y=1.1, hjust=0,label= "HA state",size=5)+
  theme_bw()+
  theme(text = element_text(size=25),
        axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0,size=20),
        plot.title = element_text(face = "bold"))+
  scale_x_datetime(date_breaks = "4 hour", labels = date_format("%H:%M"))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  ylab("Probability of every state")+ xlab("Clock time")+
  ggtitle('d') #Figure 3A.d

##Produce Figure 4A assuming the harmonic HMMs have been fitted on the physical activity time series of all the day-shift subjects
sleepAllDS  <- c() 
for (id in parameter$id[parameter$type=='Day']){
  load(paste0("HMC_",id,"_new.Rdata"))
  two_day_prob <- rbind(HMM_results_1hmc$one_day_summary_subject$one_day_prob,HMM_results_1hmc$one_day_summary_subject$one_day_prob)
  start <- which(hour(two_day_prob$time)==20 & minute(two_day_prob$time)<5)[1] #All the daily profiles start around 20:00 
  end <- start+288-1 #288=24/(1/12)
  one_day_prob <- two_day_prob$p1[start:end]
  sleepAllDS <- rbind(sleepAllDS, one_day_prob)
  rownames(sleepAllDS)[nrow(sleepAllDS)] <- paste0('S', id)
}
time <- seq(0,24-sf,sf)
sleepAllDS <- data.frame(t(sleepAllDS))
sleepAllDS <- cbind(time,sleepAllDS)
sleepAllDSforplot <- melt(sleepAllDS,id.vars = 'time', variable.name = 'series')
median_DS <- apply(sleepAllDS[,2:ncol(sleepAllDS)], 1,median)

ggplot()+
  geom_line(data=sleepAllDSforplot, aes(x=time,y=value,group = series),col="#C4961A",linetype = "dashed")+
  geom_line(aes(x=time, y=median_DS),col='black',size=1.5,linetype = "dashed")+
  scale_y_continuous(limits=c(0,1)) + 
  scale_x_discrete(breaks = seq(0,24,1), 
                   limits=c('21','','23','','1','','3','','5','','7','','9','','11','','13','','15','','17','','19','','21'))+
  theme(legend.position="none")+
  theme(text = element_text(size=20))+
  xlab("Clock Time")+
  ylab("Probability of Rest")+
  ggtitle("(a) DS")


##############################################################################################
##two-Harmonic HMMs (two circadian oscillators, one for work days and the other for free days)
##############################################################################################

##Identify work and free days for day-shift subjcts
all_date<-unique(time_date)
work_day<-unique(all_date[which(all_date %in% day(work_time))])
free_day<-all_date[-which(all_date %in% work_day)]
slot_work<-which(time_date %in% work_day) #Divide the whole study period into two periods
slot_free<-which(time_date %in% free_day)

##Identify work and free days for night-shift subjcts
Hours <- unique(format(as.POSIXct(work_time),format = "%H:%M"))
wh <- Hours[which(as.numeric(substr(Hours,1,2))>=20)] #All the night-shift work started after 20:00
wh <- unique(as.numeric(substr(wh,1,2))) 

if(length(wh)==1) { #One subject may have different work starting hours. Use the earliest one.
  work_time_start_hour<-wh[1]
  work_date<-day(work_time[which(hour(work_time)==work_time_start_hour)])
} else{
  work_time_start_hour<-wh[1]
  work_time_start_hour_2<-wh[2]
  work_date<-day(work_time[which((hour(work_time)==work_time_start_hour)|(hour(work_time)==work_time_start_hour_2)) ])
} 
work_time_start_hour <- min(wh) 

all_date<-unique(day(info$time))
free_date<-all_date[-which(all_date %in% work_date)]
work_day<-unique(work_date) #remove possible duplicates

if(hour(time[length(time)])<work_time_start_hour && free_date[length(free_date)]==time_date[length(time_date)]){
  free_day<-free_date[-length(free_date)]
} else{
  free_day<-unique(free_date)
}

if (any(diff(free_day)>1)){ #Subjects may not have consecutive free days (at most two periods of free days)
  A1 <- which(time_date %in% free_day[1])
  A2 <- which(time_date %in% free_day[which(a>1)+1])
  B <-  which(time_hour==work_time_start_hour)
  C1 <- min(which(A1%in%B))
  slot_free_start_1 <- A1[C1]
  slot_free_length_1 <- which(a>1)
  C2 <- min(which(A2%in%B))
  slot_free_start_2 <- A2[C2]
  slot_free_length_2 <- length(free_day)-which(a>1)
  slot_free <- c(slot_free_start_1:(slot_free_start_1+slot_free_length_1*24/sf-1),
                 slot_free_start_2:(slot_free_start_2+slot_free_length_2*24/sf-1))
}else{
  A<-which(time_date %in% free_day)
  B<-which(time_hour==work_time_start_hour)
  C<-min(which(A%in%B))
  slot_free_start<-A[C]
  if (free_day[length(free_day)]==all_date[length(all_date)]) {
    slot_free<-c(slot_free_start:length(time))
  } else {
    if (free_day[1]==all_date[1]){
      slot_free<-c(1:(1+length(free_day)*(24/sf)-1))
    } else{
      slot_free<-c(slot_free_start:(slot_free_start+length(free_day)*(24/sf)-1)) #index of free day
      if (slot_free[length(slot_free)]>length(time)){
        slot_free <- c(slot_free_start:length(time))
      }
    }
  }
} 
slot_work<-c(seq(1:length(time))[-slot_free]) 

##Create two circadian oscillators
lag<-0:(nrow(info)-1)
sin_baseline<-sin(2*pi*lag*sf/24)
cos_baseline<-cos(2*pi*lag*sf/24)
sin_work<-rep(0,nrow(info))
cos_work<-rep(0,nrow(info))
sin_free<-rep(0,nrow(info))
cos_free<-rep(0,nrow(info))
sin_work[slot_work]<-sin_baseline[slot_work] 
cos_work[slot_work]<-cos_baseline[slot_work]
sin_free[slot_free]<-sin_baseline[slot_free]
cos_free[slot_free]<-cos_baseline[slot_free]

##Fit 2HMM
y_2hmc<-depmix(sqValue~1,transition=~sin_work+cos_work+sin_free+cos_free,
               data=info,nstates=3,family=gaussian(),ntimes=nrow(info)) #Assume that there are 3 hidden states and that the emission density is Gaussian
AIC <- array(NA,10) #The likelihood may have local maxima, thus try 10 starting values and compare AIC
for(a in 1:10){
  set.seed(a)
  HMMfit_2hmcs<-try(fit(y_2hmc,verbose=F),TRUE)
  if(isTRUE(class(HMMfit_2hmcs)=="try-error")) 
  { next } else {AIC[a]=AIC(HMMfit_2hmcs)} 
}
AIC <- as.data.frame(AIC)
set.seed(which(AIC==min(AIC,na.rm=TRUE)))
HMMfit_2hmcs<-fit(y_2hmc,verbose=F)
HMM_results_2hmc<-twoHarmonic_HMMs(HMM=HMMfit_2hmcs,sin1=sin_work,cos1=cos_work,sin2=sin_free,cos2=cos_free) #Summarise useful parameters 

##Compute day profile and circadian parameters
circadian_states_prob<-HMM_results_2hmc$circadian_states_prob
circadian_states_prob$time<-info$time
transition_prob<-HMM_results_2hmc$transition_prob
slot_free_used<-test_24h_consecutive(slot_free,sf)
hour_day_start_range<- seq(12,20,2)
freeday<-one_day_summary(circadian_states_prob=circadian_states_prob[slot_free_used,],
                         transition_prob=transition_prob[slot_free_used,],
                         hour_day_start_range=hour_day_start_range,sf=sf, figure=F)
slot_work_used<-test_24h_consecutive(slot_work,sf)
workday<-one_day_summary(circadian_states_prob=circadian_states_prob[slot_work_used,],
                         transition_prob=transition_prob[slot_work_used,],
                         hour_day_start_range=hour_day_start_range, sf=sf,figure=F)
HMM_results_2hmc$freeday<-freeday
HMM_results_2hmc$slot_free_used<-slot_free_used
HMM_results_2hmc$workday<-workday
HMM_results_2hmc$slot_work_used<-slot_work_used
HMM_results_2hmc$sf<-sf
save(HMM_results_2hmc,file=paste0("2HMC_",id,"_new.Rdata"))


######################
###Clustering analysis
######################

##Find the number of clusters for NS and produce Figure S4 
library(factoextra) #Package for extract and visualize the output of exploratory multivariate data analyses
f_wss_NS<- fviz_nbclust(dataNS, hcut, method = "wss") #Dertemine the optimal number of clusters using the within cluster sums of squares
df_wss_NS <- as.data.frame(f_wss_NS$data)
df_wss_NS$clusters <- as.numeric(df_wss_NS$clusters)
max_wss_NS <- df_wss_NS$y[1]
ggplot(df_wss_NS)+
  scale_x_continuous(breaks=seq(0,10,1))+
  geom_point(aes(x=clusters,y=y/max_wss_NS),col='steelblue',size=2)+
  geom_line(aes(x=clusters,y=y/max_wss_NS),col='steelblue',size=1)+
  xlab("Number of Clusters")+
  ylab("% of Total Within-cluster Sum of Square")+
  ggtitle('NS')+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  theme(text = element_text(size=20)) 

##Divide NS into 3 clusters
library(TSclust) #Package for measurements of dissimilarity between time series
diss_matrix_whole_EUCL_NS<- diss(dataNS, METHOD = "EUCL") #Calculate euclidean distance between rest profiles
hclust_whole_EUCL_ward_NS <- hclust(diss_matrix_whole_EUCL_NS,method = "ward.D", members = NULL) #Agglomerative clustering with Ward¡¯s method
hclust_whole_EUCL_ward_NS <- cutree(hclust_whole_EUCL_ward_NS, k = 3) #Cut the tree at 3 clusters


#################################
###Multivariate linear regression
#################################
library('GJRM') #Package for various joint (and univariate) regression models
parameter$workRI_logodd <- log(parameter$workdayRI/(1-parameter$workdayRI))
parameter$freeRI_logodd <- log(parameter$freedayRI/(1-parameter$freedayRI))

##Determine the marginal distribution of RI (no significant difference between Gaussian and Weibull, use Gaussion)
EnvStats::qqPlot(parameter$workdayRI, distribution = "norm",estimate.params=T,add.line = TRUE)
EnvStats::qqPlot(parameter$workdayRI, distribution = "gamma",estimate.params=T,add.line = TRUE)
EnvStats::qqPlot(parameter$workdayRI, distribution = "weibull",estimate.params=T,add.line = TRUE)

EnvStats::qqPlot(parameter$freedayRI, distribution = "norm",estimate.params=T,add.line = TRUE)
EnvStats::qqPlot(parameter$freedayRI, distribution = "gamma",estimate.params=T,add.line = TRUE)
EnvStats::qqPlot(parameter$freedayRI, distribution = "weibull",estimate.params=T,add.line = TRUE)

##Select the copula function (F is the best)
out_NN_N <- gjrm(list(workRI_logodd ~ type + Age + scoreHO + N_of_years_NS+ 
                        type:Age + type:scoreHO + type:N_of_years_NS + Age:scoreHO + Age:N_of_years_NS + scoreHO:N_of_years_NS, 
                      freeRI_logodd ~ type + Age + scoreHO + N_of_years_NS+ 
                        type:Age + type:scoreHO + type:N_of_years_NS + Age:scoreHO + Age:N_of_years_NS + scoreHO:N_of_years_NS), 
                 data = parameter, margins = c("N", "N"),Model = "B",BivD='N')
out_NN_F <- gjrm(list(workRI_logodd ~ type + Age + scoreHO + N_of_years_NS+ 
                        type:Age + type:scoreHO + type:N_of_years_NS + Age:scoreHO + Age:N_of_years_NS + scoreHO:N_of_years_NS, 
                      freeRI_logodd ~ type + Age + scoreHO + N_of_years_NS+ 
                        type:Age + type:scoreHO + type:N_of_years_NS + Age:scoreHO + Age:N_of_years_NS + scoreHO:N_of_years_NS), 
                 data = parameter, margins = c("N", "N"),Model = "B",BivD='F')
out_NN_AMH <- gjrm(list(workRI_logodd ~ type + Age + scoreHO + N_of_years_NS+ 
                          type:Age + type:scoreHO + type:N_of_years_NS + Age:scoreHO + Age:N_of_years_NS + scoreHO:N_of_years_NS, 
                        freeRI_logodd ~ type + Age + scoreHO + N_of_years_NS+ 
                          type:Age + type:scoreHO + type:N_of_years_NS + Age:scoreHO + Age:N_of_years_NS + scoreHO:N_of_years_NS), 
                   data = parameter, margins = c("N", "N"),Model = "B",BivD='AMH')
out_NN_FGM <- gjrm(list(workRI_logodd ~ type + Age + scoreHO + N_of_years_NS+ 
                          type:Age + type:scoreHO + type:N_of_years_NS + Age:scoreHO + Age:N_of_years_NS + scoreHO:N_of_years_NS, 
                        freeRI_logodd ~ type + Age + scoreHO + N_of_years_NS+ 
                          type:Age + type:scoreHO + type:N_of_years_NS + Age:scoreHO + Age:N_of_years_NS + scoreHO:N_of_years_NS), 
                   data = parameter, margins = c("N", "N"),Model = "B",BivD='FGM')

AIC(out_NN_N, out_NN_F, out_NN_AMH, out_NN_FGM) #F has the lowest AIC
BIC(out_NN_N, out_NN_F, out_NN_AMH, out_NN_FGM)

##Model selected by AICc
out_NN_F <- gjrm(list(workRI_logodd ~ type + Age + scoreHO + N_of_years_NS+ 
                           type:Age + Age:scoreHO + Age:N_of_years_NS, 
                         freeRI_logodd ~ type + N_of_years_NS), 
                    data = parameter, margins = c("N","N"),Model = "B",BivD='F')

##Outliers
th <- 4/140
idx_rm <- which(cooks.distance(out_NN_F$gam1)>th | cooks.distance(out_NN_F$gam2)>th)

##Model selected by AICc without outliers
out_NN_F <- gjrm(list(workRI_logodd ~ type + Age + scoreHO + N_of_years_NS+ 
                           type:Age + type:scoreHO + Age:scoreHO + Age:N_of_years_NS, 
                         freeRI_logodd ~ N_of_years_NS), 
                    data = parameter[-idx_rm,], margins = c("N","N"),Model = "B",BivD='F')

##Produce Figure S6A
equationDS_age = function(x){
  out_NN_F$coefficients[1] + 
    (out_NN_F$coefficients[3] + out_NN_F$coefficients[8] * mean(parameter[-idx_rm,]$scoreHO) + out_NN_F$coefficients[9] * mean(parameter[-idx_rm,]$N_of_years_NS)) * x + 
    out_NN_F$coefficients[4] * mean(parameter[-idx_rm,]$scoreHO) + 
    out_NN_F$coefficients[5] * mean(parameter[-idx_rm,]$N_of_years_NS)
}
equationNS_age = function(x){
  out_NN_F$coefficients[1] + out_NN_F$coefficients[2] + 
    (out_NN_F$coefficients[3] + out_NN_F$coefficients[6] + out_NN_F$coefficients[8] * mean(parameter[-idx_rm,]$scoreHO) + out_NN_F$coefficients[9] * mean(parameter[-idx_rm,]$N_of_years_NS)) * x + 
    (out_NN_F$coefficients[4] + out_NN_F$coefficients[7]) * mean(parameter[-idx_rm,]$scoreHO) + 
    out_NN_F$coefficients[5] * mean(parameter[-idx_rm,]$N_of_years_NS)
}
ggplot(parameter[-idx_rm,],aes(y=workRI_logodd,x=Age,group=type,color=type))+geom_point()+
  stat_function(fun=equationDS_age,geom="line",color='#C4961A')+
  stat_function(fun=equationNS_age,geom="line",color='steelblue')+
  ylab('Log odds of RI on work days')+
  ggtitle('A')+
  theme_bw()+
  theme(text = element_text(size=20),plot.title = element_text(face = "bold"),legend.position="none")+
  scale_color_manual(values=c('#C4961A','steelblue'),labels=c('DS','NS'),name='Shift Type')


####################
###Spectrum analysis
####################
source('Spectrum_est_once.R')
source("Spectrum_est.R")
source("Pgram_est.R")
source("Pgram_smoothing.R")
source("Kg.fun.R")
source("Bandwidth_est.R")
source("SRslide.R")
source("Bootstrap_per.R")
source("Phase_correction.R")

info <- read.csv(paste0("agg_data_1hour_subject",id,'.csv')) #Read hourly average record 
time<- as.POSIXct(info$time, format="%Y-%m-%d %H:%M:%S")
sf=1 #Hourly means of chest surface temperature
lag<-info$lag 
temp<-info$temp_processed 
temp<- approx(lag,temp,xout=lag)$y
temp <- detrend(temp, tt = 'linear', bp = c())
spd_est_temp<-Spectrum_est_once(lag=lag, obs=temp,starting_time=starting_time, 
                                sf=sf, max_harm=3, bt_size=100,CI=90) 




