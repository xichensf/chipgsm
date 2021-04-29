rm(list = ls())
graphics.off()
library(igraph)
library(prodlim)
setwd('K562/Promoters/')
source('K562/Promoters/ChIP_GSM_functions.R')

# load paths for TF sample and input read counts 
TF_path<-read.table('Sample_input_match_list.txt', header = TRUE)
Num_TF=length(TF_path$TF_name)

# load gene annotation file hg19
hg19_RefSeq<-read.table('HG19_TSS.txt', header = TRUE)
Num_gene=length(hg19_RefSeq$gene_name)

#postive and negative 10k of TSS, 500 bps per window
segment_num=41 

# load read counts from multiple TFs at each 500 bps window at gene promoter regions, this step is slow 
multi_TFs<-loading_multi_TF_read_counts(TF_path, hg19_RefSeq)

# searching for candidate clusters, select candidate TFs, this step is slow
print('Searching for candidate clusters!')
Sampling_input_data_frame<-Candidate_cluster_searching(multi_TFs, Num_gene, segment_num, Num_TF)
print(paste('Number of candidate clusters:', toString(nrow(Sampling_input_data_frame$C_candidate_matrix))))
print(paste('Number of candidate TFs:', toString(length(Sampling_input_data_frame$Candidate_TFs))))


# distribution parameters estimation
print('Distribution fitting for each TF!')
alpha_I=matrix(0, nrow=length(Sampling_input_data_frame$Candidate_TFs), 1)
beta_I=matrix(0, nrow=length(Sampling_input_data_frame$Candidate_TFs), 1)
theta=matrix(0, nrow=length(Sampling_input_data_frame$Candidate_TFs), 1)

Xmin=10
Xupper=250

Pdfx_PowerLaw=matrix(0, nrow=length(Sampling_input_data_frame$Candidate_TFs), ncol=Xupper-Xmin+1)
Pdfi_Gamma=matrix(0, nrow=length(Sampling_input_data_frame$Candidate_TFs), ncol=Xupper)

for (t in 1:length(Sampling_input_data_frame$Candidate_TFs)){
  print(paste('Distribution fitting for:',Sampling_input_data_frame$Candidate_TFs[t]))
  # Gamma distribution fitting for background signals
  mean_Gamma<-mean(Sampling_input_data_frame$Sampling_window_Input_tags[,t])
  variance_Gamma<-var(Sampling_input_data_frame$Sampling_window_Input_tags[,t])
  alpha_I[t]=max(c(mean_Gamma^2/variance_Gamma, 1))
  beta_I[t]=min(c(variance_Gamma/mean_Gamma+1, 8))
  for (x in 1:Xupper){
    Pdfi_Gamma[t,x]=dgamma(x, shape= alpha_I[t], scale = beta_I[t])
  }
  Pdfi_Gamma[t,]=Pdfi_Gamma[t,]/sum(Pdfi_Gamma[t,])
  # Power Law distribution fitting for background signals
  est_powerlaw<-power.law.fit(Sampling_input_data_frame$Sampling_window_Sample_tags[,t], impelementation = "plfit")
  theta[t]=min(c(est_powerlaw$alpha, 3))
  for (x in Xmin:Xupper){
    Pdfx_PowerLaw[t,x-Xmin+1]=(x/Xmin)^(-theta[t])
  }
  Pdfx_PowerLaw[t,]=Pdfx_PowerLaw[t,]/sum(Pdfx_PowerLaw[t,])
  
  #plot(Pdfx_PowerLaw[t,c(1:41)]/Pdfi_Gamma[t,c(10:50)])
}

# optional, filter out some background regions
window_index=which(rowSums(Sampling_input_data_frame$Sampling_window_flag)>=2)
Sampling_input_data_frame$Sampling_window_ID = Sampling_input_data_frame$Sampling_window_ID[window_index]
Sampling_input_data_frame$Sampling_window_index = Sampling_input_data_frame$Sampling_window_index[window_index]
Sampling_input_data_frame$Sampling_window_Sample_tags = Sampling_input_data_frame$Sampling_window_Sample_tags[window_index,]
Sampling_input_data_frame$Sampling_window_Input_tags = Sampling_input_data_frame$Sampling_window_Input_tags[window_index,]
Sampling_input_data_frame$Sampling_window_FD = Sampling_input_data_frame$Sampling_window_FD[window_index,]
Sampling_input_data_frame$Sampling_window_flag = Sampling_input_data_frame$Sampling_window_flag[window_index,]
Sampling_input_data_frame$Cluster_mapping = Sampling_input_data_frame$Cluster_mapping[window_index,]
Sampling_input_data_frame$TF_count_summary =colSums(Sampling_input_data_frame$Sampling_window_Sample_tags)


print('ChIP-GSM Model parameter initialization!')
# ChIP-GSM model initialization
Gibbs_round=1000
Fold_change_factor=10

Num_TFs=length(Sampling_input_data_frame$Candidate_TFs)
Num_clusters=nrow(Sampling_input_data_frame$C_candidate_matrix)
Num_windows = length(Sampling_input_data_frame$Sampling_window_ID )

TX=matrix(0, nrow=Num_TFs, ncol=1)
TI=matrix(0, nrow=Num_TFs, ncol=1)
Binding_pattern=matrix(0, nrow=Num_windows, ncol=Num_TFs)
Cluster_pattern=matrix(0, nrow=Num_windows, ncol=Num_clusters)

# randomly assign each window or region a regulatory module
for (k in 1:Num_windows){
  candidate_cluster=which(Sampling_input_data_frame$Cluster_mapping[k,]>0)
  rand_index=sample(length(candidate_cluster), 1)
  Binding_pattern[k,]=Sampling_input_data_frame$C_candidate_matrix[candidate_cluster[rand_index[1]],]   
  Cluster_pattern[k,candidate_cluster[rand_index[1]]]=1
}


# determine the total number of read tags to be assigned to forgournd and background regions respectively for each TF
for (t in 1:Num_TFs){
  binding_index=which(Binding_pattern[,t]>0)
  Background_sampling_weight = rgamma(Num_windows, shape=alpha_I[t], scale=beta_I[t])
  Forground_sampling_weight = Background_sampling_weight
  Forground_sampling_weight[binding_index] = Fold_change_factor*Forground_sampling_weight[binding_index]
  Forground_sampling_weight = Forground_sampling_weight/(sum(Forground_sampling_weight))
  TX[t]=floor(sum(Forground_sampling_weight[binding_index])*Sampling_input_data_frame$TF_count_summary[t])
  TI[t]=Sampling_input_data_frame$TF_count_summary[t]-TX[t]    
}

# noise variance
alpha_N=2
beta_N=5
Nsigma=1/rgamma(1, shape=alpha_N, scale=1/beta_N)

Summition_binding=Binding_pattern
Summition_cluster=Cluster_pattern
Binding_pattern_old=Binding_pattern

# binding piror of distance
window_size=500
promoter_size=10000
Background_prior=matrix(1, nrow=nrow(Binding_pattern), ncol=ncol(Binding_pattern))
alpha_lambda=1
beta_lambda=1
Binding_prior<-sampling_lambda_distance(Sampling_input_data_frame, Binding_pattern_old, alpha_lambda, beta_lambda,
                                        Num_windows, Num_TFs, window_size, promoter_size)

print('Sampling starts!')
# for loop for sampling 
for (g in 1:1000){
  print(paste('Sampling round:',toString(g)))
  #***********************sampling X based on binding
  X <- sampling_X_cluster_distance(Sampling_input_data_frame, Binding_pattern_old, TX, Xmin, Nsigma, Xupper,
                                   Pdfx_PowerLaw, Binding_prior, Num_windows, Num_TFs)
  
  #***********************sampling I based on none binding
  I <- sampling_I_cluster_distance(Sampling_input_data_frame, Binding_pattern_old, TI, Nsigma, Xupper,
                                   Pdfi_Gamma, Background_prior, Num_windows, Num_TFs) 
  
  #***********************sampling noise based on sampled X and I
  Nsigma <- sampling_Nsigma_cluster(Sampling_input_data_frame, X, I, alpha_N, beta_N, Num_windows, Num_TFs)  

  # sampling clusters and update bindings
  Sampling_temp_results <- sampling_binding_cluster_distance(Sampling_input_data_frame, Binding_pattern_old,
                                                             X, Xmin, Xupper, I, Nsigma, 
                                                             Pdfx_PowerLaw, Pdfi_Gamma, 
                                                             Binding_prior, Background_prior,
                                                             Num_windows, Num_TFs, Num_clusters)
  # sampling distance parameter lambda
  Binding_prior<-sampling_lambda_distance(Sampling_input_data_frame, Sampling_temp_results$Binding_pattern,
                                          alpha_lambda, beta_lambda,
                                          Num_windows, Num_TFs, window_size, promoter_size)
  
  # control all tags used for sampling 
  Read_count_total<-sampling_TX_TI_cluster(Sampling_input_data_frame, Fold_change_factor, alpha_I, beta_I,
                                           Sampling_temp_results$Binding_pattern, Num_TFs) 
  
  TX<-Read_count_total$TX
  TI<-Read_count_total$TI
  
  # samples accumulation
  print('Propotion of all candidate positions!')
  print(sum(Sampling_temp_results$Binding_pattern)/sum(Sampling_input_data_frame$Sampling_window_flag))
  print('Propotion of previous round of candidate positions!')
  print(sum(Sampling_temp_results$Binding_pattern*Binding_pattern_old)/sum(Sampling_temp_results$Binding_pattern))
  
  Binding_pattern_old=Sampling_temp_results$Binding_pattern
  
  Summition_binding=Summition_binding+Sampling_temp_results$Binding_pattern
  Summition_cluster=Summition_cluster+Sampling_temp_results$Cluster_pattern
  if (g%%10==1){
    save(Summition_binding, Summition_cluster, file='K562_temp_results.Rdata')
  }
}


