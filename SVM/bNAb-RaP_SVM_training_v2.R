#==================================================
#==================================================
#==================================================
# Libraries
library( bio3d )
library(data.table)
library(kernlab)
library(lattice)
library(ROCR)
library(gplots)
library(caret)



#==================================================
#==================================================
#==================================================
# Main

#==================================================
# start time measurement
start.main <- proc.time()
#==================================================


#==================================================
# Command line arguments
file.training.alignment <- commandArgs()[ 3 ]
file.training.neutralization <- commandArgs()[ 4 ]
sigma2 <- as.numeric( commandArgs()[ 5 ] )
l <- as.numeric( commandArgs()[ 6 ] )
file.run.kernel <- commandArgs()[ 7 ]

# file.run.kernel <- "$PATH/run_kernelsvm.R"

#==================================================
# Source function
source( file.run.kernel )

#==================================================
# Set working directory

#==================================================
# Set up variables

#==================================================
#==================================================
# Load data

#==================================================
# Load alignment
aln <- read.fasta( file.training.alignment )
aln.ali <- aln$ali
aln.id <- aln$id

# Translate alignment into vector of strings
vec.seqs <- NULL
for( i in 1:nrow( aln.ali ) )
{
  vec.seqs <- c( vec.seqs,
                 paste( aln.ali[ i, ],
                        collapse = "" ) )
}

#==================================================
# Load neutralization data
vec.neut <- scan( file.training.neutralization )


#==================================================
#==================================================
# SVM training

#==================================================
# Get oligo kernel function and build kernel matrix
oligo_kernel1 <- oligodot( sigma2 = sigma2,
                           l = l )
kernel_matrix <- kernelMatrix( oligo_kernel1,
                               vec.seqs )
rownames( kernel_matrix ) <- aln.id
colnames( kernel_matrix ) <- aln.id

# Save kernel matrix
write.csv( kernel_matrix,
           file = "KM_training.csv",
           row.names = TRUE,
           # col.names = TRUE,
           quote = FALSE )

# Build the model using the kernel matrix with C-SVM
svm.model <- ksvm( x = kernel_matrix,
                   y = vec.neut,
                   type = "C-svc",
                   prob.model = T,
                   C = 0.01 )
save( svm.model,
      file = "svm_model.Rdata" )

##################################################
##################################################
##################################################
# For Raghav
# Determine cutoff and save it out!!!
#Get train predictions and optimal threshold using method we discussed
y_train_predict <- predict( svm.model,
                            kernel_matrix,
                            type = "p" )
pred_obj <- prediction( predictions = y_train_predict[ , 2 ],
                        labels = vec.neut )
perf_obj <- performance( pred_obj,
                         measure = "tpr", "fpr" )
pdf( "Threshold.pdf" )
plot( perf_obj,
      avg = "threshold",
      colorize = TRUE,
      lwd = 3,
      print.cutoffs.at = seq( 0, 1, by = 0.05 ),
      text.adj = c( -0.5, 0.5 ),
      text.cex = 0.6 )
grid( col = "lightgray" )
axis(1, at=seq(0, 1, by=0.1))
axis(2, at=seq(0, 1, by=0.1))
abline(v=c(0.1, 0.3, 0.5, 0.7, 0.9), col="lightgray", lty="dotted")
abline(h=c(0.1, 0.3, 0.5, 0.7, 0.9), col="lightgray", lty="dotted")
lines(x=c(0, 1), y=c(0, 1), col="black", lty="dotted")
dev.off()
#Select cutoff as minimum threshold for which FPR rate is less than 0.10
fdr_cutoff <- 0.10
fpr_ids <- which(perf_obj@x.values[[1]]<fdr_cutoff)
fpr_id <- fpr_ids[length(fpr_ids)]
opt_cutoff <- perf_obj@alpha.values[[1]][fpr_id] 
write(opt_cutoff,"Optimal_cutoff.txt")

##################################################
##################################################
##################################################





#==================================================
# Measure and print out time needed for the script
end.main <- proc.time()
duration.main <- end.main-start.main
print( paste( "Script duration:", round( duration.main[3] / 60, 2 ), "min") )
#==================================================

# Main end
#==================================================
#==================================================
#==================================================