#These functions are provided by Huang (2018).

#' oneHarmonic_HMMs
#' for shift workers that we consider two types of circadian rhythms
#' Given the output of depmix (), this function summarises some useful results for further analysis and plotting

#' @param HMM:  circadian harmonic HMM fitting results  (the return of depmix function). 3 activity states are assumed.
#' @param sin1/cos1 : one harmonics series
#' @param sf: sampling frequency

# in our model, 3 states are assumed. The order is:
# 1: inactive state; 2: moderately active state; 3: highly active state

#' @return ML_states: maximal likelihood states at each time point (results of local decoding)
#' @return prob_ML_states: probability of ML_states at each time point (results of local decoding)
#' @return transition_prob (oscillates with 24-h period):  transition probability 
#' @return circadian_states_prob  (oscillates with 24-h period):   state probability
#' @return AIC and BIC
#' @return obs_density_sq: mean and standard  deviation of  observation densities (square root), conditioned on 3 states
#' @return obs_density: 5%, 50%, 95% quantile of observation densities (original scale), conditioned on 3 states


oneHarmonic_HMMs<-function(HMM,sin1,cos1,sf=1/12){
  
  obs_params<-data.frame(mean=summary(HMM)[1:3],sd=summary(HMM)[4:6])
  
  state1_lable<-which.min(obs_params$mean)
  state3_lable<-which.max(obs_params$mean)
  state2_lable<-which(obs_params$mean==median(obs_params$mean))
  
  #################
  ###local decoding
  #################
  e1<-forwardbackward(HMM)$gamma
  L<-nrow(e1)
  state1<-e1[,state1_lable]
  state2<-e1[,state2_lable]
  state3<-e1[,state3_lable]
  states<-rep(NA,L)
  
  for(i in 1:L){
    states[i]<-ifelse((state1[i]>state2[i] & state1[i]>state3[i]),0,ifelse(state2[i]>state3[i],1,2))
  }
  
  #ML_states: maximal likelihood of current states 
  ML_states<-states+1
  #ML_states has 3 numbers where 1: inactive; 2: moderately active; 3: highly active
  
  #probability of each states
  prob_ML_states<-cbind(state1,state2,state3)
  
  #####################################################################
  ###computue the time-varing transition probabilities: transition_prob
  #####################################################################
  trans1_1<-rep(NA,L) # transition from state 1 to 1
  trans1_2<-rep(NA,L) # transition from state 1 to 2
  trans1_3<-rep(NA,L)
  trans2_1<-rep(NA,L)
  trans2_2<-rep(NA,L)
  trans2_3<-rep(NA,L)
  trans3_1<-rep(NA,L)
  trans3_2<-rep(NA,L)
  trans3_3<-rep(NA,L)
  
  A<-getpars(HMM)
  
  tmat_all<-array(NA,c(3,3,3))
  
  for (i in 1:3){
    start<-(i-1)*3+4
    end<-start+2
    tmat_all[1,i,]<-A[start:end]
  }
  for (i in 1:3){
    start<-(i-1)*3+13
    end<-start+2
    tmat_all[2,i,]<-A[start:end]
  }
  for (i in 1:3){
    start<-(i-1)*3+22
    end<-start+2
    tmat_all[3,i,]<-A[start:end]
  }
  tmat1<-tmat_all[state1_lable,,]
  tmat2<-tmat_all[state2_lable,,]
  tmat3<-tmat_all[state3_lable,,]
  
  transition_from_one_state<-function(a11,a22,a33){
    a123<-c(a11,a22,a33)
    a1<-a123[state1_lable]
    a2<-a123[state2_lable]
    a3<-a123[state3_lable]
    trans<-c(a1,a2,a3)/sum(c(a1,a2,a3))
    #inf may occour
    trans1<-replace(trans, is.na(trans), (1-sum(trans[!is.na(trans)]))/ sum(is.na(trans)) )
    return(trans1)  
  }
  
  for (i in 1:L){
    sin_1_i<-sin1[i]
    cos_1_i<-cos1[i]
    
    # transition from state 1
    a11<-exp(c(1,sin_1_i,cos_1_i)%*% tmat1[,1])
    a22<-exp(c(1,sin_1_i,cos_1_i)%*% tmat1[,2])
    a33<-exp(c(1,sin_1_i,cos_1_i)%*% tmat1[,3])
    
    y1<-transition_from_one_state(a11,a22,a33)
    trans1_1[i]<-y1[1]
    trans1_2[i]<-y1[2]
    trans1_3[i]<-y1[3]
    
    # transition from state 2
    
    a11<-exp(c(1,sin_1_i,cos_1_i)%*% tmat2[,1])
    a22<-exp(c(1,sin_1_i,cos_1_i)%*% tmat2[,2])
    a33<-exp(c(1,sin_1_i,cos_1_i)%*% tmat2[,3])
    
    y2<-transition_from_one_state(a11,a22,a33)
    trans2_1[i]<-y2[1]
    trans2_2[i]<-y2[2]
    trans2_3[i]<-y2[3]
    
    # transition from state 3
    
    a11<-exp(c(1,sin_1_i,cos_1_i)%*% tmat3[,1])
    a22<-exp(c(1,sin_1_i,cos_1_i)%*% tmat3[,2])
    a33<-exp(c(1,sin_1_i,cos_1_i)%*% tmat3[,3])
    
    y3<-transition_from_one_state(a11,a22,a33)
    trans3_1[i]<-y3[1]
    trans3_2[i]<-y3[2]
    trans3_3[i]<-y3[3]
  }
  transition_prob<-data.frame(trans1_1=trans1_1,trans1_2=trans1_2,trans1_3=trans1_3,
                              trans2_1=trans2_1,trans2_2=trans2_2,trans2_3=trans2_3,
                              trans3_1=trans3_1,trans3_2=trans3_2,trans3_3=trans3_3)
  
  
  #####################################################################
  ###computue the probability of each states (time-varing): prob_states
  #####################################################################
  P_1<-rep(NA,L)
  P_2<-rep(NA,L)
  P_3<-rep(NA,L)
  P_sum<-rep(NA,L)
  
  P_1[1]<-P_2[1]<-P_3[1]<-0
  
  P_sum[1]<-1
  initial_state<-A[1:3]
  states_0_guess<-which.max(initial_state)
  
  if(state1_lable==states_0_guess){states_0<-1;P_1[1]<-1}
  if(state2_lable==states_0_guess){states_0<-2;P_2[1]<-1}
  if(state3_lable==states_0_guess){states_0<-3;P_3[1]<-1}
  
  for (i in 2:L){
    trans_matrix<-matrix(NA,nrow=3,ncol=3)
    trans_matrix[1,]<-c(trans1_1[i],trans1_2[i],trans1_3[i])
    trans_matrix[2,]<-c(trans2_1[i],trans2_2[i],trans2_3[i])
    trans_matrix[3,]<-c(trans3_1[i],trans3_2[i],trans3_3[i])
    
    
    P_pervious<-c(P_1[i-1],P_2[i-1],P_3[i-1])
    P_1[i]<-P_pervious %*% trans_matrix[,1]
    P_2[i]<-P_pervious %*% trans_matrix[,2]
    P_3[i]<-P_pervious %*% trans_matrix[,3]
    P_sum[i]<-P_1[i]+P_2[i]+P_3[i]
  }
  
  circadian_states_prob<-data.frame(state_1=P_1,state_2=P_2,state_3=P_3)
  
  ###############observation densities#################
  mean_obs<-obs_params$mean
  
  mean_sq<-c(mean_obs[state1_lable],mean_obs[state2_lable],mean_obs[state3_lable])
  sd_sq<-c(obs_params$sd[state1_lable],obs_params$sd[state2_lable],obs_params$sd[state3_lable])
  nc_prams<-(mean_sq/sd_sq)^2
  
  var<-sd_sq^2
  obs_density<-matrix(NA,3,3)
  for (i in 1:3){
    obs_density[i,]<- qchisq(p=c(0.05,0.5,0.95), df=1, ncp = nc_prams[i])*var[i]
  }
  
  obs_density_sq<-data.frame(mean=mean_sq,sd=sd_sq)
  
  return(list(ML_states=ML_states, prob_ML_states=prob_ML_states,
              transition_prob=transition_prob,
              circadian_states_prob=circadian_states_prob,
              AIC=AIC(HMM),BIC=BIC(HMM),n_states=3,
              obs_density=obs_density,obs_density_sq=obs_density_sq))
  
}

#For twoHarmonic_HMMs.R, readers can add one more harmonic oscillator sin2 & cos2 to oneHarmonic_HMMs.R, and modify the codes according to the parameters in the 2HHMMs.


#' oneday_summary
#' day profiles summary

#' @param probs_slot: states probabilites in used slot
#' @param trans_slot : transition probabilites in used slot
#' @param hour_day_start_defined: given hour_day_start
#' @param hour_day_start_range:  a range of hour_day_start and the algorithm will choose the one with highest RI if hour_day_start is not given,
#' @param sf: sampling frequency

# in our model, 3 states are assumed. The order is:
# 1: inactive state; 2: moderately active state; 3: highly active state

#' @return one_day_prob: 3 states probabilities of one day
#' @return clocktime: corresponding time
#' @return trans_probs_oneday: transition probabilities of one day
#' @return rest_amount, center_time (gravity center of sleep), RI (rythmn index), 
# notes: when there are two phases of sleep, center_time and RI are not reliable


one_day_summary<-function(circadian_states_prob,transition_prob,
                          hour_day_start_defined=NA,hour_day_start_range=seq(12,18,2),
                          sf=1/12,figure=F){
  
  if(ncol(transition_prob)==4){  probs_states<-data.frame(time=circadian_states_prob$time,p1=circadian_states_prob$state_1,
                                                          p2=circadian_states_prob$state_2) }
  
  
  if(ncol(transition_prob)==9){
    probs_states<-data.frame(time=circadian_states_prob$time,p1=circadian_states_prob$state_1,
                             p2=circadian_states_prob$state_2,
                             p3=circadian_states_prob$state_3)  }
  
  # make periodic probs
  start<-1
  oneday_index<-start:(start+24/sf-1)
  
  probs_states_oneday<-probs_states[oneday_index,]
  probs_states_twoday<-rbind(probs_states_oneday,probs_states_oneday)
  transition_prob_oneday<-transition_prob[1:(24/sf),]
  transition_prob_twoday<-rbind(transition_prob_oneday,transition_prob_oneday)
  
  hour_time<-hour(probs_states_twoday$time)
  min_time<-minute(probs_states_twoday$time)
  
  
  find_day_start<-function(hour_day_start,hour_time,min_time){
    A<-which(hour_time==hour_day_start)
    B<-which(min_time<60*sf) 
    one_day_start<-A[min(which(A%in%B, arr.ind = TRUE))]  
    return(one_day_start)
  }
  
  select_hour_day_start<-function(hour_day_start_range){
    
    RI_range<-rep(NA,length(hour_day_start_range))
    L<-24/sf
    index<-seq(0,24-sf,sf)#absolute index to compute the gravity centre of p1
    
    for (i in 1: length(hour_day_start_range)){
      
      hour_day_start<-hour_day_start_range[i]
      one_day_start<-find_day_start(hour_day_start=hour_day_start,hour_time=hour_time,min_time=min_time)
      one_day_end<-one_day_start+24/sf-1
      
      one_day_prob<-probs_states_twoday[one_day_start:one_day_end,]
      trans_probs_oneday<-transition_prob_twoday[one_day_start:one_day_end,]
      
      index_clocktime<-seq(hour_day_start,hour_day_start+24-sf,sf)
      clocktime<-index_clocktime
      clocktime[index_clocktime>24]<-(index_clocktime-24)[index_clocktime>24]
      
      
      one_day_prob$clocktime<-clocktime
      
      # three stata probs
      p1<-one_day_prob$p1
      # duration of sleep per day
      rest_amount<-24*sum(p1)/nrow(one_day_prob)
      
      #center of rest which corresponds to the gravity centre of p1
      center_rest<-sum(p1*index/sum(p1))
      
      center_rest_clock<-clocktime[which.min(abs(index-center_rest))]
      worst_p1<-rep(rest_amount/24,24/sf)
      
      find_perfect_p1<-function(center_rest,rest_amount,index,sf){
        perfect_p1<-rep(0,24/sf)
        t1<-center_rest-rest_amount/2
        t2<-center_rest+rest_amount/2
        perfect_t2<-which.min(abs(index-center_rest))
        perfect_t1<-which.min(abs(index-center_rest+rest_amount/2))
        perfect_t3<-which.min(abs(index-center_rest-rest_amount/2))
        perfect_p1[perfect_t1:perfect_t3]<-1
        if(t2>max(index)){offset<-round((t2-max(index))/sf);perfect_p1[1:offset]<-1}
        if(t1<min(index)){offset<-round((min(index)-t1)/sf);perfect_p1[-offset+L,L]<-1}
        return(perfect_p1)
      }
      
      perfect_p1<-find_perfect_p1(center_rest,rest_amount,index,sf)
      RI_range[i]<-(sum(p1[which(perfect_p1>0)])*sf/rest_amount-rest_amount/24)*24/(24-rest_amount)
    }
    hour_day_start<-hour_day_start_range[which.is.max(RI_range)]
    return(hour_day_start)
  }
  
  if(  !is.na(hour_day_start_defined)){
    hour_day_start<-hour_day_start_defined
  }else{
    hour_day_start<-select_hour_day_start(hour_day_start_range)  }
  
  one_day_start<-find_day_start(hour_day_start=hour_day_start,hour_time=hour_time,min_time=min_time)
  one_day_end<-one_day_start+24/sf-1
  
  one_day_prob<-probs_states_twoday[one_day_start:one_day_end,]
  trans_probs_oneday<-transition_prob_twoday[one_day_start:one_day_end,]
  
  index_clocktime<-seq(hour_day_start,hour_day_start+24-sf,sf)
  clocktime<-index_clocktime
  clocktime[index_clocktime>24]<-(index_clocktime-24)[index_clocktime>24]
  
  one_day_prob$clocktime<-clocktime
  L<-24/sf
  # three stata probs
  p1<-one_day_prob$p1
  # duration of sleep per day
  rest_amount<-24*sum(p1)/nrow(one_day_prob)
  
  #center of rest which corresponds to the gravity centre of p1
  index<-seq(0,(24-sf),sf)#absolute index to compute the gravity centre of p1
  center_rest<-sum(p1*index/sum(p1))
  
  center_rest_clock<-clocktime[which.min(abs(index-center_rest))]
  worst_p1<-rep(rest_amount/24,24/sf)
  
  
  find_perfect_p1<-function(center_rest,rest_amount,index,sf){
    perfect_p1<-rep(0,24/sf)
    t1<-center_rest-rest_amount/2
    t2<-center_rest+rest_amount/2
    perfect_t2<-which.min(abs(index-center_rest))
    perfect_t1<-which.min(abs(index-center_rest+rest_amount/2))
    perfect_t3<-which.min(abs(index-center_rest-rest_amount/2))
    perfect_p1[perfect_t1:perfect_t3]<-1
    if(t2>max(index)){offset<-round((t2-max(index))/sf);perfect_p1[1:offset]<-1}
    if(t1<min(index)){offset<-round((min(index)-t1)/sf);perfect_p1[-offset+L,L]<-1}
    return(perfect_p1)
  }
  
  perfect_p1<-find_perfect_p1(center_rest,rest_amount,index,sf)
  RI<-(sum(p1[which(perfect_p1>0)])*sf/rest_amount-rest_amount/24)*24/(24-rest_amount)
  one_day_prob$perfect_p1<-perfect_p1
  one_day_prob$worst_p1<-worst_p1

  RI<-round(RI,digits=3)
  
  if(figure){
    par(mfrow=c(1,1))
    plot(index,p1,main=paste(id,'; hour start',hour_day_start,'; RI',RI))
    abline(v=center_rest,col='red',lwd=2)
  }
  clocktime[which(clocktime==24)]<-0
  p11<-trans_probs_oneday$trans1_1
  weight<-p1/sum(p1)
  weighted_p11<-sum(p11 %*% weight)
  
  one_day_summary<-list(one_day_prob=one_day_prob,clocktime=clocktime,
                        trans_probs_oneday=trans_probs_oneday,
                        rest_amount=rest_amount,weighted_p11=weighted_p11,
                        center_time=center_rest_clock,RI=RI)
  
  return(one_day_summary)
}



test_24h_consecutive<-function(slot_considered,sf){
  
  A<-diff(slot_considered)
  B<-which(A>1)
  num.slot<-length(B)+1
  slot_break_points<-matrix(NA,num.slot,3)
  
  if(length(B)>0){
    slot_break_points[1,1]<-slot_considered[1]
    slot_break_points[1,2]<-slot_considered[B[1]]
    if(length(B)>1){
      for (i in 2:length(B)){
        slot_break_points[i,1]<-slot_considered[B[i-1]+1]
        slot_break_points[i,2]<-slot_considered[B[i]]}
      slot_break_points[num.slot,1]<- slot_considered[B[i]+1]
      slot_break_points[num.slot,2]<- slot_considered[length(slot_considered)]
    }else{
      slot_break_points[num.slot,1]<- slot_considered[B[1]+1]
      slot_break_points[num.slot,2]<- slot_considered[length(slot_considered)]}
  }else{
    slot_break_points[1,1]<-slot_considered[1]
    slot_break_points[1,2]<-slot_considered[length(slot_considered)]
  }
  
  slot_break_points[,3]<-slot_break_points[,2]-slot_break_points[,1]+1
  oneday_num<-24/sf
  
  lag_points<-data.frame(lag_start=slot_break_points[,1],lag_end=slot_break_points[,2],slothours=slot_break_points[,3]*sf)
  print(lag_points)
  
  if(max(slot_break_points[,3])<oneday_num){print('no consecutive one day can be used!')
    slot_used<-NA}
  
  if(max(slot_break_points[,3])>=oneday_num){
    index1<-which.max(slot_break_points[,3])
    slot_used<-slot_break_points[index1,1]:slot_break_points[index1,2]
    slot_used<-slot_used[(length(slot_used)+1-24/sf):length(slot_used)]  #only need 24h 
  }
  
  return(slot_used)
  
}

