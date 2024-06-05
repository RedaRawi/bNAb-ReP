#==================================================
#==================================================
#==================================================
# Libraries
library( data.table )
library( h2o )
library( foreach )
library( doParallel )

library( ROCR )

library( geosphere )


#==================================================
#==================================================
#==================================================
# Main

#==================================================
# start time measurement
start.main <- proc.time()
#==================================================


#==================================================
# Set working directory

#==================================================
# Command line arguments
n.cores <- as.numeric( commandArgs()[ 3 ] )
registerDoParallel( n.cores )
memory <- as.numeric( commandArgs()[ 4 ] )
path <- commandArgs()[ 5 ]
file.training.data <- commandArgs()[ 6 ]


#==================================================
# Source function


#==================================================
# Set up variables
grid.id <- "grid"
model.id <- "final"
model.cutoff.id <- "final_cutoff"
retraining.id <- "retrain"


#==================================================
# Set up h2o
localH2O <- h2o.init( max_mem_size = paste( memory,
                                            "g",
                                            sep = "" ),
                      nthreads = n.cores,
                      enable_assertions = FALSE )

#==================================================
# Load training data
MC.path <- file.path( file.training.data )
data.hex <- h2o.importFile( path = MC.path,
                            destination_frame = "MC.hex" )
colnames( data.hex )[ length( colnames( data.hex ) ) ] <- "response"
train.names <- colnames( data.hex )[ 1:( length( colnames( data.hex ) ) - 1 ) ]

# Set response column as factor
data.hex$response <- as.factor( data.hex$response )

#==================================================
# Set GBM parameters to test
ntrees_opt <- list( 1000 )
max_depth_opt <- list( 1,
                       2,
                       3,
                       4,
                       5,
                       6 )
learnrate_opt <- list( 0.001,
                       0.01,
                       0.05,
                       0.1,
                       0.2 )
col_sample_rate_opt <- list( sqrt( length( train.names ) ) / length( train.names ),
                             0.1,
                             0.2,
                             0.3 )


hyper_parameters <- list( ntrees = ntrees_opt,
                          max_depth = max_depth_opt,
                          learn_rate = learnrate_opt,
                          col_sample_rate = col_sample_rate_opt )

#==================================================
#==================================================
# Train GBM model
gbm.model.grid <- h2o.grid( "gbm",
                            hyper_params = hyper_parameters,
                            x = train.names,
                            y = "response",
                            training_frame = data.hex,
                            nfolds = 10,
                            grid_id = grid.id,          
                            balance_classes = TRUE,
                            max_after_balance_size = 5,
                            fold_assignment = "Stratified",
                            stopping_metric = "AUC",
                            stopping_rounds = 3,
                            stopping_tolerance = 0.001 )


#==================================================
# You have to use the same metric for fetching and sorting the models.
# Fetch grid models 
grid <- h2o.getGrid( grid.id,
                     sort_by = "auc",
                     decreasing=TRUE )

#==================================================
# Sort by auc
h2o.getGrid( grid.id,
             sort_by = "auc",
             decreasing = TRUE )
print( grid )
# Save grid
save( grid, file = "Grid.Rdata")

#==================================================
# Find the best model and its full set of parameters
best_model <- h2o.getModel( grid@model_ids[[ 1 ]] )

# Get parameters
var.ntrees <- as.numeric( grid@summary_table[ 1, ]$ntrees )
var.max_depth <- as.numeric( grid@summary_table[ 1, ]$max_depth )
var.learn_rate <- as.numeric( grid@summary_table[ 1, ]$learn_rate )
var.col_sample_rate <- as.numeric( grid@summary_table[ 1, ]$col_sample_rate )

write( c( var.ntrees,
          var.max_depth,
          var.learn_rate,
          var.col_sample_rate  ),
       file = "best_parameters.txt",
       ncolumns = 1 )


#==================================================
#==================================================
# Retrain model using full data
gbm.model.final <- h2o.gbm( x = train.names,
                            y = "response",
                            training_frame = data.hex,
                            model_id = model.id,
                            balance_classes = TRUE,
                            max_after_balance_size = 5,
                            stopping_metric = "AUC",
                            stopping_rounds = 3,
                            stopping_tolerance = 0.001,
                            ntrees = var.ntrees,
                            max_depth = var.max_depth,
                            learn_rate = var.learn_rate,
                            col_sample_rate = var.col_sample_rate )

#==================================================
# Save final retrained best model 
h2o.saveModel( object = gbm.model.final,
               path = path, 
               force = TRUE )

#==================================================
# Save feature importance
df.feature.importance.raw <- as.data.frame( h2o.varimp( gbm.model.final ) )
df.feature.importance <- data.frame( cbind( df.feature.importance.raw$variable,
                                            df.feature.importance.raw$percentage ) )
colnames( df.feature.importance ) <- c( "feature",
                                        "importance" )
file.final.model.variable.importance <- paste( model.id,
                                               "_featureImportance.csv",
                                               sep = "" )
write.csv( df.feature.importance,
           file = file.final.model.variable.importance,
           row.names = FALSE,
           quote = FALSE )


#==================================================
#==================================================
# Re-train model for cutoff selection (75% / 25% split)

n.train <- round( 0.75 * nrow( data.hex ) )
n.test <- nrow( data.hex ) - n.train

vec.indices <- 1:nrow( data.hex )
vec.train <- sort( sample( vec.indices,
                           n.train ) )
# data.hex.train <- data.hex[ vec.train, ]

#==================================================
# Special case, previous one does not work with R4.3.1
df.data.hex <- as.data.frame( data.hex )
df.data.hex.train <- df.data.hex[ vec.train, ]
data.hex.train <- as.h2o( df.data.hex.train )
#==================================================

vec.test <- vec.indices[ !( vec.indices %in% vec.train ) ]
# data.hex.test <- data.hex[ vec.test, ]
#==================================================
# Special case, previous one does not work with R4.3.1
df.data.hex <- as.data.frame( data.hex )
df.data.hex.test <- df.data.hex[ vec.test, ]
data.hex.test <- as.h2o( df.data.hex.test )
#==================================================


y <- as.numeric( as.vector( data.hex[ , ncol( data.hex ) ] ) )
y.train <- y[ vec.train ]
sum( y.train == 0 ) / length( y.train )
y.test <- y[ vec.test ]
sum( y.test == 0 ) / length( y.test )

# Save training and test
h2o.exportFile( data.hex.train, 
                path = paste( path,
                              "retrain_hex_train.csv",
                              sep = "/" ) )
h2o.exportFile( data.hex.test, 
                path = paste( path,
                              "retrain_hex_test.csv",
                              sep = "/" ) )

#==================================================
# Retrain model using 75% of the data
gbm.model.final.75 <- h2o.gbm( x = train.names,
                               y = "response",
                               training_frame = data.hex.train,
                               model_id = model.cutoff.id,
                               balance_classes = TRUE,
                               max_after_balance_size = 5,
                               stopping_metric = "AUC",
                               stopping_rounds = 3,
                               stopping_tolerance = 0.001,
                               ntrees = var.ntrees,
                               max_depth = var.max_depth,
                               learn_rate = var.learn_rate,
                               col_sample_rate = var.col_sample_rate )

#==================================================
# Save final retrained best model 
h2o.saveModel( object = gbm.model.final.75,
               path = path, 
               force = TRUE )

#==================================================
# Predict using remaining 25%
y.predict.raw <- h2o.predict( gbm.model.final.75,
                              data.hex.test )
vec.raw <- as.numeric( as.vector( y.predict.raw$p1 ) )

df.predictions.true <- data.frame( cbind( vec.raw,
                                          y.test ) )
colnames( df.predictions.true ) <- c( "probabilities",
                                      "true" )

# Make predictions
mod.prediction <- prediction( df.predictions.true$probabilities,
                              df.predictions.true$true )

# Make performance
mod.performance <- performance( mod.prediction,
                                "tpr",
                                "tnr" )

#==================================================
# Extract cutoff values
vec.tnr.cutoff <- unlist( mod.performance@x.values )
vec.tpr.cutoff <- unlist( mod.performance@y.values )
vec.cutoff <- unlist( mod.performance@alpha.values )

# Set up values which you define as thresholds
cutoff.tnr <- 0.9
vec.dist <- abs( vec.tnr.cutoff - cutoff.tnr )
vec.ind.cutoff.tnr <- tail( which( vec.dist == min( vec.dist ) ), 1 )
final.cutoff <- vec.cutoff[ vec.ind.cutoff.tnr ]

#==================================================
# Calculate closest point to line
mat.bisectrix <- cbind( seq( 0, 1, 0.001 ),
                        seq( 0, 1, 0.001 ) )
mat.points <- cbind( vec.tnr.cutoff,
                     vec.tpr.cutoff )
mat.dist.bisectrix <- dist2Line( mat.points,
                                 mat.bisectrix )
mat.dist.bisectrix.sorted <- mat.dist.bisectrix[ order( mat.dist.bisectrix[ , 1 ] ), ]

df.tnr.tpr.cutoff <- data.frame( cbind( vec.tnr.cutoff,
                                        vec.tpr.cutoff,
                                        vec.cutoff ) )
colnames( df.tnr.tpr.cutoff ) <- c( "tnr",
                                    "tpr",
                                    "cutoff" )
# Select best cutoff
vec.dist.tnr <- abs( df.tnr.tpr.cutoff$tnr - as.numeric( mat.dist.bisectrix.sorted[ 1, 2 ] ) )
# df.tnr.tpr.cutoff$tnr[ which( vec.dist.tnr == min( vec.dist.tnr ) ) ]
vec.dist.tpr <- abs( df.tnr.tpr.cutoff$tpr - as.numeric( mat.dist.bisectrix.sorted[ 1, 3 ] ) )
# df.tnr.tpr.cutoff$tpr[ which( vec.dist.tpr == min( vec.dist.tpr ) ) ]

final.cutoff.balanced <- df.tnr.tpr.cutoff[ df.tnr.tpr.cutoff$tnr == unique( df.tnr.tpr.cutoff$tnr[ which( vec.dist.tnr == min( vec.dist.tnr ) ) ] ) &
                                              df.tnr.tpr.cutoff$tpr == unique( df.tnr.tpr.cutoff$tpr[ which( vec.dist.tpr == min( vec.dist.tpr ) ) ] ), ]$cutoff

#==================================================
# Save final cutoff probability
write( final.cutoff,
       file = paste( model.cutoff.id,
                     "_cutoff.txt",
                     sep = "" ) )
write( final.cutoff.balanced,
       file = paste( model.cutoff.id,
                     "_cutoff_balanced.txt",
                     sep = "" ) )

#==================================================
# TNR vs. TPR plot
pdf( "TNRvsTPR.pdf",
     width = 8,
     height = 8 )
par( mar = c( 5, 
              5,
              4,
              2 ) )

plot( mod.performance,
      avg = "threshold",
      colorize = TRUE,
      lwd = 3,
      print.cutoffs.at = seq( 0, 
                              1,
                              by = 0.05 ),
      text.adj = c( -0.5,
                    0.5 ),
      text.cex = 0.6 )
grid( col = "lightgray" )
axis( 1, 
      at = seq( 0,
                1, 
                by = 0.1 ) )
axis( 2,
      at = seq( 0,
                1,
                by = 0.1 ) )
abline( v = c( 0.1,
               0.3,
               0.5,
               0.7,
               0.9 ),
        col = "lightgray",
        lty = "dotted" )
abline( h = c( 0.1,
               0.3,
               0.5,
               0.7,
               0.9 ),
        col = "lightgray",
        lty = "dotted" )
lines( x = c( 0,
              1 ),
       y = c( 0,
              1 ),
       col = "black",
       lty = "dotted" )
dev.off()




#==================================================
#==================================================
# Retrain 10 times
n <- 10
for( i in 1:n )
{
  print( "==================================================" )
  print( paste( i,
                n,
                sep = "/" ) )
  
  retraining.id.i <- paste( retraining.id,
                            "_",
                            i,
                            sep = "" )
  
  gbm.model.retrain <- h2o.gbm( x = train.names,
                                y = "response",
                                training_frame = data.hex,
                                model_id = retraining.id.i,
                                balance_classes = TRUE,
                                nfolds = 10,
                                max_after_balance_size = 5,
                                stopping_metric = "AUC",
                                stopping_rounds = 3,
                                stopping_tolerance = 0.001,
                                ntrees = var.ntrees,
                                max_depth = var.max_depth,
                                learn_rate = var.learn_rate,
                                col_sample_rate = var.col_sample_rate )
  df.performance.i <- as.data.frame( gbm.model.retrain@model$cross_validation_metrics_summary )
  
  # Save performance
  write.csv( df.performance.i,
             file = paste( retraining.id.i,
                           ".csv",
                           sep = "" ),
             quote = FALSE )
}


#==================================================
# Close h2o
h2o.shutdown( prompt = FALSE )





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