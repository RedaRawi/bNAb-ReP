#==================================================
#==================================================
#==================================================
# Libraries
library( bio3d )
library( data.table )
library( kernlab )
library( lattice )
library( ROCR )
library( gplots )
library( caret )
library( MLmetrics )



#==================================================
#==================================================
#==================================================
# Function


#==================================================
# Convert probabilities
get.probabilities <- function( vec.raw,
                               threshold )
{
  # get_class_one_predictions <- as.numeric(as.vector(y.predict.raw$p1))
  
  # Convert the raw score into a probability
  y.predict.prob <- rep( 0,
                         length( vec.raw ) )
  for (i in 1:length( vec.raw ) )
  {
    if (vec.raw[ i ] >= threshold )
    {
      y.predict.prob[ i ] <- 0.5 + 0.5 * ( ( vec.raw[ i ] - threshold ) / ( 1 - threshold ) )
    }
    else
    {
      y.predict.prob[ i ] <- 0.5 - 0.5 * ( ( threshold - vec.raw[ i ] ) / threshold )
    }
  }
  return( y.predict.prob )
}



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
sigma2 <- as.numeric( commandArgs()[ 3 ] )
l <- as.numeric( commandArgs()[ 4 ] )
file.KM.training <- commandArgs()[ 5 ]
file.alignment.training <- commandArgs()[ 6 ]
file.neut.training <- commandArgs()[ 7 ]
file.alignment.testing <- commandArgs()[ 8 ]
file.neut.testing <- commandArgs()[ 9 ]
file.run.kernel <- commandArgs()[ 10 ]


#==================================================
# Source function
source( file.run.kernel )
load( file.training.model )

#==================================================
# Set working directory

#==================================================
# Set up variables
file.optimal.cutoff <- "Optimal_cutoff.txt"
opt_cutoff <- scan( file.optimal.cutoff )


#==================================================
#==================================================
# Load data

#==================================================
# i) Load training alignment
aln.training <- read.fasta( file.training.alignment )
aln.training.ali <- aln.training$ali
aln.training.id <- aln.training$id
vec.seqs.training <- NULL
for( i in 1:nrow( aln.training.ali ) )
{
  vec.seqs.training <- c( vec.seqs.training,
                          paste( aln.training.ali[ i, ],
                                 collapse = "" ) )
}

#==================================================
# Align testing sequences to training alignment
cmd.align <- paste( "mafft",
                    "--add",
                    file.testing.alignment,
                    "--keeplength",
                    file.training.alignment,
                    ">",
                    "Alignment.fasta" )
system( cmd.align )

#==================================================
# Extract testing sequences

# Load aligned sequences
aln <- read.fasta( "Alignment.fasta" )
aln.ali <- aln$ali
aln.id <- aln$id

# Extract aligned test sequences
aln.ali.testing <- aln.ali[ (nrow( aln.training.ali ) + 1):nrow( aln.ali ), ]
aln.id.testing <- aln.id[ (nrow( aln.training.ali ) + 1):nrow( aln.ali ) ]
vec.seqs.testing <- NULL
for( i in 1:nrow( aln.ali.testing ) )
  # for( i in 1:x )
{
  vec.seqs.testing <- c( vec.seqs.testing,
                         paste( aln.ali.testing[ i, ],
                                collapse = "" ) )
}

#==================================================
# Load neutralization data
y.test <- scan( file.testing.neutralization )


#==================================================
#==================================================
# SVM model prediction

##################################################
##################################################
##################################################
# Doubl-check with RM if this correct???
oligo_kernel1 <- oligodot( sigma2 = sigma2,
                           l = l )
kernel_matrix <- kernelMatrix( oligo_kernel1,
                               vec.seqs.testing,
                               vec.seqs.training )
rownames( kernel_matrix ) <- aln.id.testing
colnames( kernel_matrix ) <- aln.training.id
# Save kernel matrix
write.csv( kernel_matrix,
           file = "KM_testing.csv",
           row.names = TRUE,
           # col.names = TRUE,
           quote = FALSE )


# Get probabilities 
vec.probabilities <- predict( svm.model,
                              kernel_matrix,
                              type = "probabilities" )[ , 2 ]

# Save probabilities
write( vec.probabilities,
       file = "SVM-probabilities.txt",
       ncolumns = 1 )

# ##################################################
# # @RR: Should we not transform the probabilities using this function get.probabilities()???
# vec.probabilities <- get.probabilities( vec.probabilities,
#                                         opt_cutoff )
# y.predict <- rep( 0,
#                  length( vec.probabilities ) )
# y.predict[ vec.probabilities > 0.5 ] <- 1
# ##################################################

# Get opt_cutoff from Optimal threshold file
y.predict <- ifelse( vec.probabilities >= opt_cutoff, 1, 0 )
pred_test_obj <- prediction( predictions = y.predict, 
                             labels = y.test )
perf_test_obj <- performance( pred_test_obj,
                              measure = "auc" )

# Get auc for each fold (this is equivalent to balanced accuracy from confusionmatrix of caret)
auc <- round( perf_test_obj@y.values[[ 1 ]],
              4 )
print( paste( "AUC:",
              round( auc,
                     4 ) ) )

##################################################
##################################################
##################################################


#==================================================
#==================================================
# Access prediction accuracy (e.g. AUC, ACC, MSE, RMSE)
df.predictions <- data.frame( cbind( y.predict,
                                     vec.probabilities ) )
colnames( df.predictions ) <- c( "predict",
                                 "p1" )
df.predictions.true <- data.frame( cbind( df.predictions,
                                          y.test ) )
colnames( df.predictions.true ) <- c( colnames( df.predictions ),
                                      "true" )
print( df.predictions.true )

#==================================================
# AUC
# if( sum( df.predictions.true$true == 0 ) > 0 )
# {
#   mod.prediction <- prediction( df.predictions.true$p1,
#                                 df.predictions.true$true )
#   mod.auc <- performance( mod.prediction,
#                           measure = "auc" )
#   auc <- round( as.numeric( mod.auc@y.values ), 4 )
#   print( paste( "auc: ",
#                 auc ) )
# }

#==================================================
# acc
acc <- round( sum( df.predictions.true$true == df.predictions.true$predict ) / length( df.predictions.true$true ),
              4 )
print( paste( "Accuracy:",
              round( acc,
                     4 ) ) )

#==================================================
# mcc
mcc <- get_mcc( S = df.predictions.true$predict,
                Y = df.predictions.true$true )
print( paste( "MCC:",
              round( mcc,
                     4 ) ) )

#==================================================
# rmse
rmse <- round( RMSE( df.predictions.true$predict,
                     df.predictions.true$true ),
               4 )
print( paste( "rmse: ",
              rmse ) )
rmse <- round( RMSE( df.predictions.true$p1,
                     df.predictions.true$true ),
               4 )
print( paste( "rmse: ",
              rmse ) )

#==================================================
# mse
mse <- round( MSE( df.predictions.true$predict,
                   df.predictions.true$true ),
              4 )
print( paste( "mse: ",
              mse ) )

if( sum( df.predictions.true$true == 0 ) > 0 )
{
  df.metrics <- data.frame( cbind( c( "AUC",
                                      "ACC",
                                      "RMSE",
                                      "MSE" ),
                                   c( auc,
                                      acc,
                                      rmse,
                                      mse ) ) )
} else
{
  df.metrics <- data.frame( cbind( c( "ACC",
                                      "RMSE",
                                      "MSE" ),
                                   c( acc,
                                      rmse,
                                      mse ) ) )
}
colnames( df.metrics ) <- c( "metrics",
                             "value" )
print( df.metrics )

# write.csv( df.metrics,
#            # file = "FinalPerformance_SVM.csv",
#            file = "FinalPerformance_SVM-Hake.csv",
#            row.names = FALSE,
#            quote = FALSE )





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