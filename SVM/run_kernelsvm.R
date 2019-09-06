library(data.table)
library(kernlab)
library(Matrix)
library(proxy)
library(pracma)
warnings('off')

#Find if a substring exists in a string
mfunc <- function(sub_string,seq)
{
  return(unlist(gregexpr(pattern=sub_string,text=seq)))
}

#Find all substrings in order of length n in a string x
allsubstr <- function(x, n) 
{
  return(unique(substring(x, 1:(nchar(x) - n + 1), n:nchar(x))))
}

#Calculate exponential of large negative numbers
myexp <- function(y)
{
  return(ifelse(y<=-40,0,exp(y)))
}

#Oligo kernel implementation
calculate_olig_kernel <- function(seq1,seq2,subseq_l,var_param)
{
  all_subseq_length_l_seq1 <- allsubstr(seq1,subseq_l)
  all_subseq_length_l_seq2 <- allsubstr(seq2,subseq_l)
  intersect_subseqs <- intersect(all_subseq_length_l_seq1,all_subseq_length_l_seq2)
  p_list <- sapply(intersect_subseqs,mfunc,seq1)
  q_list <- sapply(intersect_subseqs,mfunc,seq2)
  sim_value <- 0
  for (i in 1:length(p_list))
  {
    #if (unlist(gregexpr("--",names(p_list[i])))[1]==-1)
    {
      p_info <- as.numeric(unlist(p_list[[i]]))
      q_info <- as.numeric(unlist(q_list[[i]]))
      distances <- as.vector(dist(p_info,q_info)^2)
      all_values <- (-1*distances)/(4*var_param)
      all_imp_values <- sapply(all_values,myexp)
      sim_value <- sim_value+sum(all_imp_values)
    }
  }
  sim_value <- sqrt(pi)*sqrt(var_param)*sim_value
  return(sim_value)
}

#Kernel generating function
oligodot <- function(sigma2=1,l=2){
  oligo_kernel <- function(x,y){
    if (is(x,"character") && is(y,"character") && nchar(x)==nchar(y) && strcmp(x,y)) {return(1)}
    else {
      out_xx <- calculate_olig_kernel(x,x,l,sigma2)
      out_yy <- calculate_olig_kernel(y,y,l,sigma2)
      out_xy <- calculate_olig_kernel(x,y,l,sigma2)
      kval <- out_xy/(sqrt(out_xx)*sqrt(out_yy));
      return(kval)
    }
  }
  return (new ("oligo_kernel",.Data=oligo_kernel, kpar=list(sigma,l)))
}
setClass("oligo_kernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

#oligo_kernel1 <- oligodot(3.6,2)
#oligo_kernel1(get_data$Env_Align_Seq[1],get_data$Env_Align_Seq[2])
#new_kernel_matrix <- kernelMatrix(oligo_kernel1,get_data$Env_Align_Seq)