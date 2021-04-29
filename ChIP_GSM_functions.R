#******************************************************************
#                         ChIP-GSM functions
#                             Xi Chen,
#                 CBIL, ECE Department,  Virginia Tech
#                            09/29/2015
#
#*******************************************************************

loading_multi_TF_read_counts<-function(TF_path, hg19_RefSeq){
  Num_TF=length(TF_path$TF_name)
  Num_gene=length(hg19_RefSeq$gene_name)
  segment_num=41 #postive and negative 10k of TSS, 500 bps per window
  
  Y_Region_ID=matrix(nrow=Num_gene*segment_num, ncol=1)
  Y_Segment_ID=matrix(nrow=Num_gene*segment_num, ncol=1)
  Y_Input_tags=matrix(0, nrow=Num_gene*segment_num, ncol=Num_TF)
  Y_Sample_tags=matrix(0, nrow=Num_gene*segment_num, ncol=Num_TF)
  
  Normalize_T=5e6
  Y_upper=250
  
  for (j in 1:Num_TF){
    print(sprintf('%s', TF_path$TF_name[j]))
    Region_Sample_tags_full=matrix(0, nrow=length(hg19_RefSeq$gene_name), ncol=segment_num)
    Region_Input_tags_full=matrix(0, nrow=length(hg19_RefSeq$gene_name), ncol=segment_num)
    
    TF_sample_count<-as.matrix(read.table(file=sprintf("%s",TF_path$sample_file[j]), row.names=1, header = TRUE))
    Gene_sample_list<-rownames(TF_sample_count)
    index<-match(Gene_sample_list, hg19_RefSeq$gene_name)
    for (i in 1:length(index)){
      Region_Sample_tags_full[index[i],]=as.numeric(TF_sample_count[i,]) 
    }
    
    TF_input_count<-as.matrix(read.table(file=sprintf("%s",TF_path$input_file[j]), row.names=1, header = TRUE))
    Gene_input_list<-rownames(TF_input_count)  
    index<-match(Gene_input_list, hg19_RefSeq$gene_name)
    for (i in 1:length(index)){
      Region_Input_tags_full[index[i],]=as.numeric(TF_input_count[i,])
    }
    
    for (k in 1:Num_gene){
      for (s in 1:segment_num){
        Y_Region_ID[(k-1)*segment_num+s,1]=hg19_RefSeq$gene_name[k]
        Y_Segment_ID[(k-1)*segment_num+s,1]=s-21
        Y_Input_tags[(k-1)*segment_num+s,j]=round(Region_Input_tags_full[k,s])
        Y_Sample_tags[(k-1)*segment_num+s,j]=round(Region_Sample_tags_full[k,s])
      }
    }
    
    Y_Sample_tags[which(Y_Sample_tags<0),j]=0
    Y_Input_tags[which(Y_Input_tags<0),j]=0
    
    Y_Sample_tags[,j]=floor(Normalize_T/(sum(Y_Sample_tags[,j]))*Y_Sample_tags[,j])
    Y_Input_tags[,j]=floor(Normalize_T/(sum(Y_Input_tags[,j]))*Y_Input_tags[,j])
  }
  
  Y_Input_tags[which(Y_Input_tags>=Y_upper)]=Y_upper
  Y_Sample_tags[which(Y_Sample_tags>=Y_upper)]=Y_upper
  FD=(Y_Sample_tags+1)/(Y_Input_tags+1)  
  
  multi_TFs=list("Y_Region_ID"=Y_Region_ID, 
                 "Y_Segment_ID"=Y_Segment_ID,
                 "Y_Sample_tags"=Y_Sample_tags,
                 "Y_Input_tags"=Y_Input_tags, 
                 "FD"=FD)
  return(multi_TFs)
}
#******************************************************************







#*******************************************************************
Candidate_cluster_searching<-function(multi_TFs, Num_gene, segment_num, Num_TF){
  Sampling_window_flag=matrix(1, nrow=Num_gene*segment_num, ncol=Num_TF)
  Sampling_window_flag[which(multi_TFs$FD<2, arr.ind=T)]=0
  Sampling_window_flag[which(multi_TFs$Y_Sample_tags<50, arr.ind=T)]=0
  Sampling_window_flag=Sampling_window_flag[which(rowSums(Sampling_window_flag)<=5),]
  Sampling_window_flag=Sampling_window_flag[which(rowSums(Sampling_window_flag)>=2),]
  Cluster_pattern<-unique(Sampling_window_flag)
  Candidate_cluster_num=matrix(0, nrow=nrow(Cluster_pattern), ncol=1)
  #index=which(row.match(data.frame(Sampling_window_flag), data.frame(Cluster_pattern[1,]))>0)
  #index=which(apply(Sampling_window_flag, 1, identical, Cluster_pattern[1,]))
  index=which(rowSums(abs(sweep(Sampling_window_flag, 2, Cluster_pattern[1,])))==0)
  if (sum(index)>0){
    Candidate_cluster_num[1,1]=length(index)
  } 
  for (k in 2:nrow(Cluster_pattern)){
      #index=which(row.match(data.frame(Sampling_window_flag), data.frame(Cluster_pattern[k,]))>0)
      #index=which(apply(Sampling_window_flag, 1, identical, Cluster_pattern[k,]))
      index=which(rowSums(abs(sweep(Sampling_window_flag, 2, Cluster_pattern[k,])))==0)
      if (sum(index)>0){
        Candidate_cluster_num[k,1]=length(index)
        Sampling_window_flag[index,]=0
        Sampling_window_flag=Sampling_window_flag[which(rowSums(Sampling_window_flag)>0),]
      } 
  }

  
  cluster_index=which(Candidate_cluster_num>50)
  TF_index=which(colSums(Cluster_pattern[cluster_index,])>0)
  Candidate_TFs=TF_path$TF_name[TF_index]
  C_candidate_matrix=Cluster_pattern[cluster_index,TF_index]
  C_candidate_matrix=rbind(C_candidate_matrix, 0)
  Total_clusters=nrow(C_candidate_matrix)
  
  Sampling_window_flag=matrix(1, nrow=Num_gene*segment_num, ncol=length(TF_index))
  Sampling_window_ID=multi_TFs$Y_Region_ID
  Sampling_window_index=multi_TFs$Y_Segment_ID
  Sampling_window_Sample_tags=multi_TFs$Y_Sample_tags[,TF_index]
  Sampling_window_Input_tags=multi_TFs$Y_Input_tags[,TF_index]
  Sampling_window_FD=multi_TFs$FD[,TF_index]
  Sampling_window_flag[which(Sampling_window_FD<2)]=0
  Sampling_window_flag[which(Sampling_window_Sample_tags<10)]=0
  TF_count_summary=colSums(Sampling_window_Sample_tags)
  
  Cluster_mapping=matrix(0, nrow=Num_gene*segment_num, ncol=Total_clusters)
  for (k in 1:nrow(Sampling_window_flag)){
    Cluster_mapping[k,Total_clusters]=1
    if (sum(Sampling_window_flag[k,])>=2){
      for (s in 1:Total_clusters-1){
        if (sum(abs(C_candidate_matrix[s,]-Sampling_window_flag[k,]*C_candidate_matrix[s,]))==0){
          Cluster_mapping[k,s]=1
        }
      }
    }
  }

  Sampling_input_data_frame=list('Sampling_window_ID'=Sampling_window_ID,
                                 'Sampling_window_index'=Sampling_window_index,
                                 'Sampling_window_Sample_tags'=Sampling_window_Sample_tags,
                                 'Sampling_window_Input_tags'=Sampling_window_Input_tags,
                                 'Sampling_window_FD'=Sampling_window_FD,
                                 'Sampling_window_flag'=Sampling_window_flag,
                                 'Cluster_mapping'=Cluster_mapping,
                                 'C_candidate_matrix'=C_candidate_matrix,
                                 'Candidate_TFs'=Candidate_TFs,
                                 'TF_count_summary'=TF_count_summary)
  return(Sampling_input_data_frame)
}
#*******************************************************************






#*******************************************************************
sampling_X_cluster_distance<-function(Sampling_input_data_frame, Binding_pattern, TX, Xmin, Nsigma, Xupper, Pdfx_PowerLaw, Binding_prior, Num_windows, Num_TFs){
X=matrix(0, nrow=Num_windows, ncol=Num_TFs)
Gaussion_pdf=matrix(0, nrow=2*Xupper+1, ncol=1)
for (k in -Xupper:Xupper){
  Gaussion_pdf[k+Xupper+1,1]=dnorm(k, mean=0, sd=sqrt(Nsigma))
}
Gaussion_pdf=Gaussion_pdf/sum(Gaussion_pdf)

for (t in 1:Num_TFs){
  PX=matrix(0, nrow=Num_windows, ncol=1)
  # assign Xmin tags to each binding region
  binding_index= which(Binding_pattern[,t]>0)
  Number_bindings=length(binding_index)
  X[binding_index, t]=Xmin
  # determine the weight for each binding region
  for (k in 1:Number_bindings){
    mean1=Sampling_input_data_frame$Sampling_window_Sample_tags[binding_index[k], t]
    Y_min=max(c(floor(mean1-3*sqrt(Nsigma)), Xmin))
    Y_max=min(c(floor(mean1+3*sqrt(Nsigma)), Xupper))   
    if (Y_min<Y_max){
      Pdfx_Gaussion=Gaussion_pdf[c((Xupper+Y_min-mean1):(Xupper+Y_max-mean1))]
      Pdfx_Gaussion=Pdfx_Gaussion/sum(Pdfx_Gaussion)
      
      Pdfx=Pdfx_Gaussion*Pdfx_PowerLaw[t,c((Y_min-Xmin+1):(Y_max-Xmin+1))]
      Pdfx=Pdfx/sum(Pdfx)
      
      X_temp = sample(c(Y_min:Y_max), 1, prob=Pdfx)
    }else{
      X_temp=Y_min#Y_min=Ymax, only one choice
    }
    PX[binding_index[k], 1] = Pdfx_PowerLaw[t,X_temp-Xmin+1]*dnorm(X_temp, mean=mean1, sd=sqrt(Nsigma))*Binding_prior[binding_index[k],t]
  }
  
  # assign tags to each region according to the weight defined in PX
  PX_tag=PX[binding_index,1]
  PX_tag=PX_tag/sum(PX_tag)
  Tag_location= sample(c(1:Number_bindings), TX[t]-Xmin*Number_bindings, replace = TRUE, prob=PX_tag)
  X_add<-hist(Tag_location, breaks=c(0.5:(Number_bindings+0.5)), plot = FALSE)
  #     X_add=floor(X_add/100)
  X[binding_index, t]=X[binding_index, t]+X_add$counts 
}

X[which(X>Xupper)]=Xupper
 return(X)
}
#*******************************************************************






#*******************************************************************
sampling_I_cluster_distance<-function(Sampling_input_data_frame, Binding_pattern, TI, Nsigma, Xupper, Pdfi_Gamma, Background_prior, Num_windows, Num_TFs){
I=matrix(0, nrow=Num_windows, ncol=Num_TFs)

Gaussion_pdf=matrix(0, nrow=2*Xupper+1, ncol=1)
for (k in -Xupper:Xupper){
  Gaussion_pdf[k+Xupper+1,1]=dnorm(k, mean=0, sd=sqrt(Nsigma))
}
Gaussion_pdf=Gaussion_pdf/sum(Gaussion_pdf)

for (t in 1:Num_TFs){
  PI=matrix(0, nrow=Num_windows, ncol=1)
  # assign Xmin tags to each binding region
  None_binding_index = which(Binding_pattern[,t]==0)
  Num_none_binding = length(None_binding_index)
  # determine the weight for each binding region
  for (k in 1:Num_none_binding){
    mean1=Sampling_input_data_frame$Sampling_window_Sample_tags[None_binding_index[k], t]
    Y_min=max(c(floor(mean1-3*sqrt(Nsigma)), 1))
    Y_max=min(c(floor(mean1+3*sqrt(Nsigma)), Xupper))
    if (Y_min < Y_max){
      Pdfi_Gaussion=Gaussion_pdf[c((Xupper+Y_min-mean1):(Xupper+Y_max-mean1))]
      Pdfi_Gaussion=Pdfi_Gaussion/sum(Pdfi_Gaussion)
      Pdfi=Pdfi_Gaussion*Pdfi_Gamma[t, c(Y_min:Y_max)]
      Pdfi=Pdfi/sum(Pdfi)
      I_temp = sample(c(Y_min:Y_max),1,prob=Pdfi) 
    }else{
      I_temp=Y_min
    }
    PI[None_binding_index[k]] = Pdfi_Gamma[t, I_temp]*dnorm(I_temp, mean=mean1, sd=sqrt(Nsigma))*Background_prior[None_binding_index[k], t]
  }
  # assign tags to each region according to the weight defined in PI
  PI_tag=PI[None_binding_index, 1]
  PI_tag=PI_tag/sum(PI_tag)
  
  Tag_location= sample(c(1:Num_none_binding), TI[t], replace=TRUE, prob=PI_tag)
  I_add<-hist(Tag_location, breaks=c(0.5:(Num_none_binding+0.5)), plot = FALSE)
  I[None_binding_index, t]=I[None_binding_index, t]+I_add$counts 
} 
I[which(I>Xupper)]=Xupper
return(I)
}
#*******************************************************************






#*******************************************************************
sampling_Nsigma_cluster<-function(Sampling_input_data_frame, X, I, alpha_N, beta_N, Num_windows, Num_TFs){
  alpha_new=alpha_N+1/2
  beta_new=beta_N+sum((Sampling_input_data_frame$Sampling_window_Sample_tags-X-I)^2)/(2*Num_windows*Num_TFs)
  Nsigma = 1/rgamma(1, shape=alpha_new, scale=1/beta_new)
  #Nsigma = beta_new/rchisq(1, 2*alpha_new)
  return(Nsigma) 
}  
#********************************************************************




#********************************************************************
sampling_binding_cluster_distance<-function(Sampling_input_data_frame, Binding_pattern_old, X, Xmin, Xupper, I, Nsigma,
                                  Pdfx_PowerLaw, Pdfi_Gamma, Binding_prior, Background_prior, Num_windows, Num_TFs, Num_clusters){
  Binding_pattern=matrix(0, nrow=Num_windows, ncol=Num_TFs)
  Cluster_pattern=matrix(0, nrow=Num_windows, ncol=Num_clusters)
  Gaussion_pdf=matrix(0, nrow=2*Xupper+1, ncol=1)
  for (k in -Xupper:Xupper){
    Gaussion_pdf[k+Xupper+1,1]=dnorm(k, mean=0, sd=sqrt(Nsigma))
  }
  Gaussion_pdf=Gaussion_pdf/sum(Gaussion_pdf)
  
  for (k in 1:Num_windows){
    candiate_cluster_index=which(Sampling_input_data_frame$Cluster_mapping[k,]>0)
    if (length(candiate_cluster_index)==1){#only have zero binding, then, no sampling is needed
      #current region has no binding with sample/input > 2
      Binding_pattern[k,]=0 
      Cluster_pattern[k,Num_clusters]=1
    }
    else{
      PX_new=matrix(0, nrow=1, ncol=Num_TFs)
      PI_new=matrix(0, nrow=1, ncol=Num_TFs)
      all_possible_binding_vector=colSums(Sampling_input_data_frame$C_candidate_matrix[candiate_cluster_index,])
      all_possible_binding_vector[which(all_possible_binding_vector>0)]=1
      for (t in 1:Num_TFs){
        if (Binding_pattern_old[k,t]==0 && all_possible_binding_vector[t]==1){
          # sample X
          mean_x=Sampling_input_data_frame$Sampling_window_Sample_tags[k,t]
          Y_min=max(c(floor(mean_x-3*sqrt(Nsigma)), Xmin))
          Y_max=min(c(floor(mean_x+3*sqrt(Nsigma)), Xupper))
          if (Y_min<Y_max){
            Pdfx_Gaussion= Gaussion_pdf[c((Xupper+Y_min-mean_x):(Xupper+Y_max-mean_x))]
            Pdfx_Gaussion=Pdfx_Gaussion/sum(Pdfx_Gaussion)
            Pdfx=Pdfx_Gaussion*Pdfx_PowerLaw[t,c((Y_min-Xmin+1):(Y_max-Xmin+1))]
            Pdfx=Pdfx/sum(Pdfx)
            X_temp = sample(c(Y_min:Y_max),1,prob=Pdfx) 
          }else{
            X_temp=Y_min
          }
          PX = Pdfx_PowerLaw[t, X_temp-Xmin+1]*dnorm(X_temp, mean=mean_x, sd=sqrt(Nsigma))*Binding_prior[k,t] 
          if (I[k,t]==0){
            PI = Pdfi_Gamma[t, I[k,t]+1]*dnorm(I[k,t], mean=mean_x, sd=sqrt(Nsigma))*Background_prior[k,t] 
          }else{
            PI = Pdfi_Gamma[t, I[k,t]]*dnorm(I[k,t], mean=mean_x, sd=sqrt(Nsigma))*Background_prior[k,t] 
          }
          PX_new[1,t] = PX/(PI+PX)
          PX_new[1,t] = max(c(min(c(0.98, PX_new[1,t])), 0.02))                    
          PI_new[1,t] = 1-PX_new[1,t]
          
        }else if(Binding_pattern_old[k,t]==1 && all_possible_binding_vector[t]==0){
          
          #sample I
          mean_i=Sampling_input_data_frame$Sampling_window_Sample_tags[k,t]
          Y_min = max(c(floor(mean_i-3*sqrt(Nsigma)), 1))
          Y_max = min(c(floor(mean_i+3*sqrt(Nsigma)), Xupper)) 
          if (Y_min<Y_max){
            Pdfi_Gaussion = Gaussion_pdf[c((Xupper+Y_min-mean_i):(Xupper+Y_max-mean_i))]
            Pdfi_Gaussion = Pdfi_Gaussion/sum(Pdfi_Gaussion)
            Pdfi = Pdfi_Gaussion*Pdfi_Gamma[t, c(Y_min:Y_max)]
            Pdfi = Pdfi/(sum(Pdfi))
            I_temp = sample(c(Y_min:Y_max),1,prob=Pdfi) 
          }else{
            I_temp=Y_min
          }
          PI = Pdfi_Gamma[t, I_temp]*dnorm(I_temp, mean=mean_i, sd=sqrt(Nsigma))*Background_prior[k,t]  
          PX = Pdfx_PowerLaw[t, X[k,t]-Xmin+1]*dnorm(X[k,t], mean=mean_i, sd=sqrt(Nsigma))*Binding_prior[k,t]   
          PX_new[1,t] = PX/(PI+PX)
          PX_new[1,t] = max(c(min(c(0.98, PX_new[1,t])), 0.02))                    
          PI_new[1,t] = 1-PX_new[1,t]
        }else {
          #all_possible_binding_vector(1,t)==0
          PX_new[1,t]=0.5
          PI_new[1,t]=0.5
        }
      }
      
      P_candidate_cluster=matrix(1,nrow=1, ncol=length(candiate_cluster_index))
      
      for (c in 1:length(candiate_cluster_index)){
        for (t in 1:Num_TFs){
          P_candidate_cluster[1,c]=P_candidate_cluster[1,c]*(PX_new[1,t]^Sampling_input_data_frame$C_candidate_matrix[candiate_cluster_index[c],t])*(PI_new[1,t]^(1-Sampling_input_data_frame$C_candidate_matrix[candiate_cluster_index[c],t]))
        }
      }
      P_candidate_cluster=P_candidate_cluster/sum(P_candidate_cluster)
      
      cluster_sample_index = sample(length(candiate_cluster_index),1,prob=P_candidate_cluster)
      Binding_pattern[k,]=Sampling_input_data_frame$C_candidate_matrix[candiate_cluster_index[cluster_sample_index],]
      Cluster_pattern[k,candiate_cluster_index[cluster_sample_index]]=1
    }
  }
  return(list('Binding_pattern'=Binding_pattern, 'Cluster_pattern'=Cluster_pattern))
}
#********************************************************






#********************************************************
sampling_TX_TI_cluster<-function(Sampling_input_data_frame, Fold_change_factor, alpha_I, beta_I, Binding_pattern, Num_TFs){
  TX=matrix(0, nrow=Num_TFs, ncol=1)
  TI=matrix(0, nrow=Num_TFs, ncol=1)
  for (t in 1:Num_TFs){
    binding_index=which(Binding_pattern[,t]>0)
    Background_sampling_weight = rgamma(Num_windows, shape=alpha_I[t], scale=beta_I[t])
    Forground_sampling_weight = Background_sampling_weight
    Forground_sampling_weight[binding_index]=Fold_change_factor*Forground_sampling_weight[binding_index]
    Forground_sampling_weight=Forground_sampling_weight/(sum(Forground_sampling_weight))
    TX[t]=floor(sum(Forground_sampling_weight[binding_index])*Sampling_input_data_frame$TF_count_summary[t]);
    TI[t]=Sampling_input_data_frame$TF_count_summary[t]-TX[t];
  }   
  return(list('TX'=TX, 'TI'=TI)) 
}
#**********************************************************




#*********************************************************
sampling_lambda_distance<-function(Sampling_input_data_frame, Binding_pattern_old, alpha_lambda, beta_lambda,
                                            Num_windows, Num_TFs, window_size, promoter_size){
  Binding_prior=matrix(0, nrow=Num_windows, ncol=Num_TFs)
  alpha_new=matrix(0, nrow=Num_TFs, ncol=1)
  beta_new=matrix(0, nrow=Num_TFs, ncol=1)
  Normalization_C=matrix(0, nrow=Num_TFs, ncol=1)
  lambda_new=matrix(0, nrow=Num_TFs, ncol=1)
  weight=matrix(0, nrow=Num_TFs, ncol=41)
  for (t in 1:Num_TFs){
    alpha_new[t]=alpha_lambda+sum(Binding_pattern_old[,t])
    beta_new[t]=beta_lambda+sum(((abs(Sampling_input_data_frame$Sampling_window_index)+0.5)*Binding_pattern_old[,t]))
    #lambda_new[t]=rgamma(1, shape = alpha_new[t], scale=1/beta_new[t])
    #lambda_new[t]=beta_new[t]/rchisq(1, 2*alpha_new[t])
    lambda_new[t]=1/(sum(((abs(Sampling_input_data_frame$Sampling_window_index)+0.5)*Binding_pattern_old[,t]))/sum(Binding_pattern_old[,t]))
    Normalization_C[t]=sum(lambda_new[t]*exp(-lambda_new[t]*c(0.5:20.5)))
    Binding_prior[,t]=(promoter_size/window_size)/Normalization_C[t]*lambda_new[t]*exp(-lambda_new[t]*(abs(Sampling_input_data_frame$Sampling_window_index)+0.5))
    #for (w in 1:41){
    #  weight[t, w]=sum(Sampling_input_data_frame$Sampling_window_Sample_tags[which(Sampling_input_data_frame$Sampling_window_index==(w-21)), t])
    #}
    #weight[t,]=41*weight[t,]/sum(weight[t,])
    #for (w in 1:41){
    #  index=which(Sampling_input_data_frame$Sampling_window_index==(w-21))
    #  Binding_prior[index,t]=weight[t,w]
    #}
  }  
  return(Binding_prior) 
}
#*********************************************************








#*********************************************************
sampling_lambda_distance_constant<-function(Sampling_input_data_frame, Binding_pattern_old, alpha_lambda, beta_lambda,
                                        Num_windows, Num_TFs, window_size, promoter_size){
Binding_prior=matrix(0, nrow=Num_windows, ncol=Num_TFs)
#alpha_new=matrix(0, nrow=Num_TFs, ncol=1)
#beta_new=matrix(0, nrow=Num_TFs, ncol=1)
#Normalization_C=matrix(0, nrow=Num_TFs, ncol=1)
#lambda_new=matrix(0, nrow=Num_TFs, ncol=1)
weight=matrix(0, nrow=Num_TFs, ncol=41)
for (t in 1:Num_TFs){
  #alpha_new[t]=alpha_lambda+sum(Binding_pattern_old[,t])
  #beta_new[t]=beta_lambda+sum(((abs(Sampling_input_data_frame$Sampling_window_index)+0.5)*Binding_pattern_old[,t]))
  #lambda_new[t]=rgamma(1, shape = alpha_new[t], scale=1/beta_new[t])
  #lambda_new[t]=1/(sum(((abs(Sampling_input_data_frame$Sampling_window_index)+0.5)*Binding_pattern_old[,t]))/sum(Binding_pattern_old[,t]))
  #Normalization_C[t]=sum(lambda_new[t]*exp(-lambda_new[t]*c(0.5:20.5)))
  #Binding_prior[,t]=(promoter_size/window_size)/Normalization_C[t]*lambda_new[t]*exp(-lambda_new[t]*(abs(Sampling_input_data_frame$Sampling_window_index)+0.5))
  for (w in 1:41){
    weight[t, w]=sum(Sampling_input_data_frame$Sampling_window_Sample_tags[which(Sampling_input_data_frame$Sampling_window_index==(w-21)), t])
  }
  weight[t,]=41*weight[t,]/sum(weight[t,])
  for (w in 1:41){
    index=which(Sampling_input_data_frame$Sampling_window_index==(w-21))
    Binding_prior[index,t]=weight[t,w]
  }
}  
Binding_prior=matrix(1, nrow=Num_windows, ncol=Num_TFs)
return(Binding_prior) 
}
#***********************************************************