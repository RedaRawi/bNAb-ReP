#==================================================
#==================================================
#==================================================
# Libraries

library( bio3d )



#==================================================
#==================================================
#==================================================
# Functions

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
#==================================================
#==================================================
# Main

#==================================================
# start time measurement
start.main <- proc.time()
#==================================================

#==================================================
# Command line arguments
file.neutralization <- commandArgs()[ 3 ]
file.alignment <- commandArgs()[ 4 ]
mafft.exe.path <- commandArgs()[ 5 ]

# file.neutralization <- "VRC01_IC50_50_neutralization.txt"
# file.alignment <- "VRC01_IC50_50_alignment.fasta"
# mafft.exe.path <- "/usr/local/bin/mafft"

#==================================================
# Source function

#==================================================
# Set working directory
# setwd( "/Users/rawir/Documents/work/vrc/projects/bNAb-ReP/data/generation-1/tmp/" )


#==================================================
# Set variables
file.training <- paste( "Training_",
                        unlist( strsplit( file.alignment,
                                          "\\_alignment.fasta" ) ),
                        ".txt",
                        sep = "" )
file.alignment.hxb2 <- paste( unlist( strsplit( file.alignment,
                                                "\\.fasta" ) ),
                              "_hxb2.fasta",
                              sep = "" )
file.alignment.hxb2.glycan <- paste( unlist( strsplit( file.alignment.hxb2,
                                                       "\\.fasta" ) ),
                                     "_glycan.fasta",
                                     sep = "" )
file.alignment.glycan <- paste( unlist( strsplit( file.alignment,
                                                  "\\.fasta" ) ),
                                "_glycan.fasta",
                                sep = "" )

AA.gap <- unlist( strsplit("ACDEFGHIKLMNPQRSTVWY-",split = "" ) )
AA.gap.X <- unlist( strsplit("ACDEFGHIKLMNPQRSTVWY-X",split = "" ) )

hxb2 <- "MRVKEKYQHLWRWGWRWGTMLLGMLMICSATEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTPLCVSLKCTDLKNDTNTNSSSGRMIMEKGEIKNCSFNISTSIRGKVQKEYAFFYKLDIIPIDNDTTSYKLTSCNTSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEINCTRPNNNTRKRIRIQRGPGRAFVTIGKIGNMRQAHCNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFYCNSTQLFNSTWFNSTWSTEGSNNTEGSDTITLPCRIKQIINMWQKVGKAMYAPPISGQIRCSSNITGLLLTRDGGNSNNESEIFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTKAKRRVVQREKRAVGIGALFLGFLGAAGSTMGAASMTLTVQARQLLSGIVQQQNNLLRAIEAQQHLLQLTVWGIKQLQARILAVERYLKDQQLLGIWGCSGKLICTTAVPWNASWSNKSLEQIWNHTTWMEWDREINNYTSLIHSLIEESQNQQEKNEQELLELDKWASLWNWFNITNWLWYIKLFIMIVGGLVGLRIVFAVLSIVNRVRQGYSPLSFQTHLPTPRGPDRPEGIEEEGGERDRDRSIRLVNGSLALIWDDLRSLCLFSYHRLRDLLLIVTRIVELLGRRGWEALKYWWNLLQYWSQELKNSAVSLLNATAIAVAEGTDRVIEVVQGACRAIRHIPRRIRQGLERILL"
file.hxb2 <- "HXB2.fasta"
write.fasta( ids = "B.FR.83.HXB2_LAI_IIIB_BRU.K03455",
             seqs = unlist( strsplit( hxb2,
                                      "" ) ),
             file = file.hxb2 )


#==================================================
# Load neutralization data
vec.neut <- scan( file.neutralization )


#==================================================
# Align HXB2 sequence to reference alignment using MAFFT
command <- paste( mafft.exe.path,
                  "--addfull",
                  file.hxb2,
                  "--keeplength",
                  file.alignment,
                  ">",
                  file.alignment.hxb2 )
system( command )

# Load MAFFT aligned sequences
aln.alignment.hxb2 <- read.fasta( file.alignment.hxb2 )
aln.alignment.hxb2.ali <- aln.alignment.hxb2$ali
aln.alignment.hxb2.id <- aln.alignment.hxb2$id



#==================================================
# Check alignment - Does it contain non-standard AA
if( sum( !( as.vector( aln.alignment.hxb2.ali ) %in% AA.gap.X ) ) > 0 )
{
  stop( "Alignment does contain non-standard amino acids. Please double-check!" )
}
#==================================================

#==================================================
# Convert glycan positions
aln.alignment.hxb2.ali.glycan <- convert.N.2.NGlyc( aln.alignment.hxb2.ali )

# Save alignment WITH HXB2
write.fasta( ids = aln.alignment.hxb2.id,
             seqs = aln.alignment.hxb2.ali.glycan,
             file = file.alignment.hxb2.glycan )

# Save alignment WITHOUT HXB2
write.fasta( ids = aln.alignment.hxb2.id[ 1:( length( aln.alignment.hxb2.id ) - 1 ) ],
             seqs = aln.alignment.hxb2.ali.glycan[ 1:( nrow( aln.alignment.hxb2.ali.glycan ) - 1 ), ],
             file = file.alignment.glycan )


#==================================================
# Perform on-hot encoding AA21
print( "Perform one-hot encoding using AA21 (20 AA, glycan)" )
df.one.hot.encoding.AA21 <- one.hot.encoding.AA21( file.alignment.glycan )
#==================================================

df.training <- data.frame( cbind( df.one.hot.encoding.AA21,
                                  vec.neut ) )
# Save training data
write.table( df.training,
             file = file.training,
             row.names = FALSE,
             col.names = FALSE,
             quote = FALSE )





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