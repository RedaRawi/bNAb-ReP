#==================================================
#==================================================
#==================================================
# Libraries

library( bio3d )
library( h2o )



#==================================================
#==================================================
#==================================================
# Functions

#==================================================
# Align test sequence to reference/training alignment using MAFFT
align.test.to.train <- function( file.test.alignment,
                                 file.reference.alignment,
                                 mafft.exe.path )
{
  #==================================================
  # Set variables
  file.test.alignment.aligned.2.reference <- paste( unlist( strsplit( file.test.alignment,
                                                                      "\\.fasta" ) ),
                                                    "_aligned2Reference.fasta",
                                                    sep = "" )
  # Number of test sequences
  n <- nrow( read.fasta( file.test.alignment, rm.dup = FALSE )$ali )
  file.test.alignment.aligned <- paste( unlist( strsplit( file.test.alignment,
                                                          "\\.fasta" ) ),
                                        "_aligned.fasta",
                                        sep = "" )
  
  #==================================================
  # Align test sequence to reference alignment using MAFFT
  command <- paste( mafft.exe.path,
                    "--quiet",
                    "--addfull",
                    file.test.alignment,
                    "--keeplength",
                    file.reference.alignment,
                    ">",
                    file.test.alignment.aligned.2.reference )
  system( command )
  
  # Load MAFFT aligned sequences
  aln.test.aligned <- read.fasta( file.test.alignment.aligned.2.reference, rm.dup = FALSE )
  aln.test.aligned.ali <- aln.test.aligned$ali
  aln.test.aligned.id <- aln.test.aligned$id
  
  # Extract only Lynch clones
  aln.test.aligned.ali <- aln.test.aligned.ali[ tail( 1:nrow( aln.test.aligned.ali ),
                                                      n ), ]
  aln.test.aligned.id <- aln.test.aligned.id[ tail( 1:length( aln.test.aligned.id ),
                                                    n ) ]
  
  #==================================================
  # Save aligned test sequences
  write.fasta( seqs = aln.test.aligned.ali,
               ids = aln.test.aligned.id,
               file = file.test.alignment.aligned )
  
  #==================================================
  # Return final alignment file name
  return( file.test.alignment.aligned )
}
#==================================================

#==================================================
# Check whether test sequences are used for training
test.sequence.in.reference.alignment <- function( file.test.alignment.aligned,
                                                  file.reference.alignment )
{
  #==================================================
  # Load reference alignment and test sequences and
  # check if test sequence in reference alignment???
  aln.test <- read.fasta( file.test.alignment.aligned, rm.dup = FALSE )
  aln.test.ali <- aln.test$ali
  aln.test.id <- aln.test$id
  
  vec.seqs.test <- NULL
  for( i in 1:nrow( aln.test.ali ) )
  {
    vec.seq <- aln.test.ali[ i, ]
    var.seq <- paste( vec.seq[ vec.seq != "-" ],
                      collapse = "" )
    vec.seqs.test <- c( vec.seqs.test,
                        var.seq )
  }
  
  #==================================================
  # Load reference alignment
  aln.reference <- read.fasta( file.reference.alignment, rm.dup = FALSE )
  aln.reference.ali <- aln.reference$ali
  vec.seqs.reference <- NULL
  for( i in 1:nrow( aln.reference.ali ) )
  {
    vec.seq <- aln.reference.ali[ i, ]
    var.seq <- paste( vec.seq[ vec.seq != "-" ],
                      collapse = "" )
    vec.seqs.reference <- c( vec.seqs.reference,
                             var.seq )
  }
  
  #==================================================
  # Test if sequences are in reference alignment 
  vec.indices.test.in.reference <- which( vec.seqs.test %in% vec.seqs.reference )
  if( length( vec.indices.test.in.reference ) > 0 )
  {
    print( "Test sequences in training set!" )
  } else
  {
    print( "Test sequences not in training set!" )
  }
  
  return( vec.indices.test.in.reference )
}
#==================================================

#==================================================
# Convert regular alignment to 21 alphabet alignment (20 standard AA plus Glycan position)
convert.N.2.NGlyc <- function( alignment )
{
  alignment.glycan <- matrix( 0,
                              nrow( alignment ),
                              ncol( alignment ) )
  
  for( i in 1:nrow( alignment ) )
  {
    vec.seq.aligned <- alignment[ i, ]
    vec.seq <- vec.seq.aligned[ vec.seq.aligned != "-" ]
    df.ind.aligned.raw <- data.frame( cbind( which( vec.seq.aligned != "-" ),
                                             1:sum( vec.seq.aligned != "-" ) ) )
    colnames( df.ind.aligned.raw ) <- c( "alignment",
                                         "raw" )
    
    vec.ind.N <- which( vec.seq == "N" )
    vec.ind.help <- NULL
    for( ind in vec.ind.N )
    {
      if( vec.seq[ ind + 2  ] %in% c( "S",
                                      "T" ) )
      {
        vec.ind.help <- c( vec.ind.help,
                           1 )
      } else
      {
        vec.ind.help <- c( vec.ind.help,
                           0 )
      }
    }
    vec.ind.glycan <- vec.ind.N[ vec.ind.help == 1 ]
    
    vec.cols.glycan <- df.ind.aligned.raw[ df.ind.aligned.raw$raw %in% vec.ind.glycan, ]$alignment
    # Assign new sequence
    alignment.glycan[ i, ] <- alignment[ i, ]
    alignment.glycan[ i, vec.cols.glycan ] <- "@"
  }
  
  return( alignment.glycan )
}
#==================================================

#==================================================
# Generate alignment with glycosylation sites annotated
alignment.2.annotated.alignment <- function( file.test.alignment.aligned )
{
  #==================================================
  # Load alignment
  aln <- read.fasta( file.test.alignment.aligned, rm.dup = FALSE )
  aln.ali <- aln$ali
  aln.id <- aln$id
  
  #==================================================
  # Convert glycan positions
  aln.ali.glycan <- convert.N.2.NGlyc( aln.ali )
  
  # Save alignment WITH HXB2
  file.test.alignment.aligned.glycan <- paste( unlist( strsplit( file.test.alignment.aligned,
                                                                 "\\.fasta" ) ),
                                               "_glycan.fasta",
                                               sep = "" )
  #==================================================
  # Save alignment
  write.fasta( ids = aln.id,
               seqs = aln.ali.glycan,
               file = file.test.alignment.aligned.glycan )
  
  #==================================================
  # Retruen alignment file name
  return( file.test.alignment.aligned.glycan )
}

#==================================================
# One-hot encoding of sequences
one.hot.encoding.AA21 <- function( file.alignment,
                                   AA.glycan = unlist( strsplit("ACDEFGHIKLMNPQRSTVWY@", split = "" ) ) )
{
  #==================================================
  # Load alignment
  aln <- read.fasta( file.alignment )
  aln.ali <- aln$ali
  aln.id <- aln$id
  
  p.AA <- length( AA.glycan )
  p <- ncol( aln.ali )
  matrix.alignment.binary <- matrix( 0,
                                     nrow( aln.ali ),
                                     ncol( aln.ali ) * p.AA )
  var.window.start <- seq( 1,
                           ( ( p - 1 ) * p.AA ) + 1,
                           p.AA)
  var.window.end <- seq( p.AA,
                         p  * p.AA,
                         p.AA )
  for( i in 1:nrow( aln.ali ) )
  {
    for( j in 1:ncol( aln.ali ) )
    {
      aa <- aln.ali[ i, j ]
      vec.col.ind <- var.window.start[ j ]:var.window.end[ j ]
      var.col <- vec.col.ind[ which( aa == AA.glycan ) ]
      matrix.alignment.binary[ i, var.col ] <- 1
    }
  }
  
  # Return converted alignment
  return( matrix.alignment.binary )
}
#==================================================

#==================================================
# Predict resistance/sensitivity probabilities
bNAb.RaP.predict <- function( path.file.gbm.model,
                              path.file.testing.features,
                              cutoff )
{
  #==================================================
  # Set up h2o
  localH2O <- h2o.init()
  
  #==================================================
  # Load best model
  gbm.model <- h2o.loadModel( path = path.file.gbm.model )
  # print( gbm.model )
  
  #==================================================
  # Load test data (H2O)
  MC.path <- file.path( path.file.testing.features )
  data.hex <- h2o.importFile( path = MC.path,
                              destination_frame = "MC.hex" )
  
  #==================================================
  # Extract predictions
  y.predict.raw <- h2o.predict( gbm.model,
                                data.hex )
  vec.raw <- as.numeric( as.vector( y.predict.raw$p1 ) )
  
  #==================================================
  # Convert to adjusted probabilities (TNR 0.9)
  vec.prob <- get.probabilities( vec.raw,
                                 cutoff )
  
  #==================================================
  # Close h2o
  h2o.shutdown( prompt = FALSE )
  
  #==================================================
  # Return prediction probabilities
  return( vec.prob )
}
#==================================================

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
# Main bNAb-ReP function
bNAb.ReP <- function( file.test.alignment,
                      file.reference.alignment,
                      mafft.exe.path,
                      path.file.gbm.model,
                      cutoff,
                      bNAb = "VRC01" )
{
  #==================================================
  # Align test sequences to reference alignment
  print( "Align test sequence(s) to reference alignment" )
  file.test.alignment.aligned <- align.test.to.train( file.test.alignment,
                                                      file.reference.alignment,
                                                      mafft.exe.path )
  
  #==================================================
  # Check if aligned test sequences in reference alignment
  print( "Check if aligned test sequence(s) in reference alignment" )
  var.test.in.training <- test.sequence.in.reference.alignment( file.test.alignment.aligned,
                                                                file.reference.alignment )
  if( sum( var.test.in.training ) > 0 )
  {
    print( paste( "Test sequence(s)", which( var.test.in.training == 1 ), "used during training!" ) )
  }
  
  #==================================================
  # Generate features
  print( "Generate features" )
  
  # Generate alignment with glycan
  print( "Annotate glycosylation positions" )
  file.test.alignment.glycan <- alignment.2.annotated.alignment( file.test.alignment.aligned )
  
  # Perform on-hot encoding AA21
  print( "Perform on-hot encoding using AA21 (20 AA, glycan)" )
  df.one.hot.encoding.AA21 <- one.hot.encoding.AA21( file.test.alignment.glycan )
  
  # Save feature file
  df.testing <- data.frame( df.one.hot.encoding.AA21 )
  write.table( df.testing,
               path.file.testing.features,
               row.names = FALSE,
               col.names = FALSE,
               quote = FALSE )
  
  #==================================================
  # Predict probability for a sequence to be sensitive
  print( "Predict probability for sequence(s) to be sensitive" )
  vec.prediction.probabilities <- bNAb.RaP.predict( path.file.gbm.model,
                                                    path.file.testing.features,
                                                    cutoff )
  
  #==================================================
  # Return probabilities
  return( vec.prediction.probabilities )
}
#==================================================



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
bNAb <- commandArgs()[ 3 ]
file.test.alignment <- commandArgs()[ 4 ]
mafft.exe.path <- commandArgs()[ 5 ]
output.prefix <- commandArgs()[ 6 ]

# #==================================================
# setwd( "/Users/rawir/Documents/work/vrc/sbis_projects/bNAb-RaP/data/generation-1/full/tmp" )
# bNAb <- "VRC34.01"
# file.test.alignment <- "tmp.fasta"
# mafft.exe.path <- "/usr/local/bin/mafft"
# #==================================================

#==================================================
# Source function

#==================================================
# Set working directory

#==================================================
# Set variables
file.reference.alignment <- paste( "alignments/", bNAb, "_IC50_50_alignment.fasta", sep = "" )
path.file.gbm.model <- paste( getwd(), "/models/", bNAb, "_full", sep = "" )
AA.gap <- unlist( strsplit("ACDEFGHIKLMNPQRSTVWY-",split = "" ) )
df.bNAb.cutoff <- data.frame( cbind( c( "10-1074", "2F5", "2G12", "35O22", "3BNC117", "4E10", "8ANC195", "CH01", "DH270.1", "DH270.5", "DH270.6", "HJ16", "NIH45-46", "PG16", "PG9", "PGDM1400", "PGT121", "PGT128", "PGT135", "PGT145", "PGT151", "VRC-CH31", "VRC-PG04", "VRC01", "VRC03", "VRC07", "VRC13", "VRC26.08", "VRC26.25", "VRC29.03", "VRC34.01", "VRC38.01", "b12" ),
                                     c( 0.6599304, 0.3763821, 0.2487729, 0.5453792, 0.9692675, 0.9689797, 0.7105667, 0.559615, 0.5539247, 0.7228051, 0.6883996, 0.3709901, 0.8690571, 0.8229699, 0.772257, 0.9187894, 0.8666024, 0.6032068, 0.2800005, 0.9009995, 0.8522265, 0.9371531, 0.8612042, 0.9408964, 0.5167682, 0.9283123, 0.7873914, 0.5849593, 0.6087503, 0.2214183, 0.5212775, 0.3016291, 0.2591896 ) ) )
colnames( df.bNAb.cutoff ) <- c( "bNAb", "cutoff" )
df.bNAb.cutoff$bNAb <- as.character( df.bNAb.cutoff$bNAb )
df.bNAb.cutoff$cutoff <- as.numeric( as.character( df.bNAb.cutoff$cutoff ) )
cutoff <- df.bNAb.cutoff[ df.bNAb.cutoff$bNAb == bNAb, ]$cutoff
if( is.na( output.prefix ) )
{
  file.testing <- paste( "Testing_",
                         as.numeric( Sys.time() ), "_",
                         bNAb, 
                         ".txt",
                         sep = "" )
  path.file.testing.features <- paste( getwd(),
                                       file.testing,
                                       sep = "/" )
  file.prediction.probabilities <- paste( "tmp_", as.numeric( Sys.time() ), "_probabilities.txt", sep = "" )
  
} else
{
  file.testing <- paste( "Testing_",
                         output.prefix, "_",
                         bNAb, 
                         ".txt",
                         sep = "" )
  path.file.testing.features <- paste( getwd(),
                                       file.testing,
                                       sep = "/" )
  file.prediction.probabilities <- paste( output.prefix, "_probabilities.txt", sep = "" )
}




#==================================================
#==================================================
#==================================================
# Start bNAb-ReP function
vec.prediction.probabilities <- bNAb.ReP( file.test.alignment,
                                          file.reference.alignment,
                                          mafft.exe.path,
                                          path.file.gbm.model,
                                          cutoff,
                                          bNAb )

# Save probabilities
write( vec.prediction.probabilities,
       file = file.prediction.probabilities,
       ncolumns = 1 )
#==================================================
#==================================================
#==================================================





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