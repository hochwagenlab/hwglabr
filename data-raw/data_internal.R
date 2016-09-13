#------------------------------------------------------------------------------#
#                                                                              #
#                 Generation of data internal to the package                   #
#                                                                              #
#------------------------------------------------------------------------------#

# Get all data frames and then generate internal package data
# at the end of the script (/R/sysdata.rda)


#------------------------------------------------------------------------------#
#             First iteration of convergent intergenic regions                 #
#                            (Now deprecated)                                  #
# S288C genome
# Import data file, originally generated using:
# "/Volumes/LabShare/Luis/LabWork/Scripts/2015.10_H3k79me3_transcription.R"
path <- '/Users/luis/Google_Drive_NYU/LabShare_Luis/LabWork/2015.05_dot1_rec8dot1/2015.10_H3K79me3_transcripts/'
S288C_conv <- read.table(paste0(path, 'conv_midpoints_rpkm.txt'),
                         header = TRUE, stringsAsFactors = FALSE)
# Remove info about transcription:
S288C_conv <- S288C_conv[, 1:6]

# SK1 genome
# Import data file, originally generated using:
# /Volumes/LabShare/Luis/LabWork/Scripts/2015.11_convMidpointsSK1.R
path <- '/Users/luis/Google_Drive_NYU/LabShare_Luis/LabWork/GenomeSequences/SK1/'
SK1_conv <- read.table(paste0(path, 'conv_midpoints_SK1.txt'),
                       header = TRUE, stringsAsFactors = FALSE)


#------------------------------------------------------------------------------#
#                       New intergenic region data                             #
#                         (Conv, div and tandem)                               #
# The intergenic region coordinate files were generated (by Luis and Tovah) using
# the scripts found at: '/Volumes/LabShare/Luis/LabWork/GenomeSequences/hwglabr/'
# 1. Import data file
path <- '/Volumes/LabShare/Luis/LabWork/GenomeSequences/hwglabr/'
S288C_conv_midpoint_dist <- read.table(paste0(path, 'S288C_conv_midpoint_dist.txt'),
                                       header = TRUE, stringsAsFactors = FALSE)
S288C_div_midpoint_dist <- read.table(paste0(path, 'S288C_div_midpoint_dist.txt'),
                                      header = TRUE, stringsAsFactors = FALSE)
S288C_tand_midpoint_dist <- read.table(paste0(path, 'S288C_tand_midpoint_dist.txt'),
                                       header = TRUE, stringsAsFactors = FALSE)
SK1_conv_midpoint_dist <- read.table(paste0(path, 'SK1_conv_midpoint_dist.txt'),
                                     header = TRUE, stringsAsFactors = FALSE)
SK1_div_midpoint_dist <- read.table(paste0(path, 'SK1_div_midpoint_dist.txt'),
                                    header = TRUE, stringsAsFactors = FALSE)
SK1_tand_midpoint_dist <- read.table(paste0(path, 'SK1_tand_midpoint_dist.txt'),
                                     header = TRUE, stringsAsFactors = FALSE)


#------------------------------------------------------------------------------#
#                               Centromeres                                    #
# SK1 info based on Keeney lab genome sequence and annotation
SK1cen <- data.frame("Chromosome" = c("chr01","chr02","chr03","chr04","chr05",
                                     "chr06","chr07","chr08","chr09","chr10",
                                     "chr11","chr12","chr13","chr14","chr15",
                                     "chr16"),
                     "Start" = c(137832, 226711, 128699, 463204, 157003, 162815,
                                 505440, 95031, 346215, 415648, 452723, 137738,
                                 249103, 616840, 307236, 553355),
                     "End" = c(137948, 226826, 128779, 463321, 157119, 162931,
                               505558, 95147, 346330, 415764, 452838, 137855,
                               249221, 616956, 307353, 553467),
                     "Mid" = c(137890, 226768, 128739, 463262, 157061, 162873,
                               505499, 95089, 346272, 415706, 452780, 137796,
                               249162, 616898, 307294, 553411),
                     "LenChr" = c(203893, 794508, 342718, 1490682, 602514,
                                  284456, 1067526, 544538, 435585, 719294,
                                  687260, 1008248, 908607, 812465, 1054033,
                                  921188))

# S288C
path <- '/Users/luis/Google_Drive_NYU/LabShare_Luis/LabWork/GenomeSequences/'
S288C_gff <- read.table(paste0(path, 'saccharomyces_cerevisiae_R64-1-1_20110208.gff'),
                        fill = TRUE, stringsAsFactors = FALSE)
S288C_gff <- S288C_gff[1:16406, ]
S288Ccen <- S288C_gff[S288C_gff[, 3] == 'centromere', c(1, 4:5)]
names(S288Ccen) <- c('Chromosome', 'Start', 'End')
S288Ccen$Mid <- floor(S288Ccen$Start +  (S288Ccen$End - S288Ccen$Start) / 2)
S288Ccen$LenChr <- S288C_gff[S288C_gff[, 3] == 'chromosome', 5][1:16]


#------------------------------------------------------------------------------#
#                               SKI rosetta                                    #
path <- '/Users/luis/Google_Drive_NYU/LabShare_Luis/LabWork/GenomeSequences/'
SK1rosetta <- read.table(paste0(path, 'SK1rosetta_R64.txt'),
                         header = TRUE, sep = '\t')

#------------------------------------------------------------------------------#
#                                GFF files                                     #
#               Used by opening_act() to run aignal_at_orf()                   #
s288C_gff <- '/Volumes/LabShare/GenomeSequences/s288C_annotation_R64_modified.gff'
SK1_gff <- '/Volumes/LabShare/GenomeSequences/SK1_MvO_V1___GENOME/SK1_annotation/SK1_annotation_modified_v2.gff'
s288C_gff <- hwglabr::gff_read(s288C_gff)
# Further parse attributes field to get gene ID only
SK1_gff <- hwglabr::gff_read(SK1_gff)
SK1_gff$attributes <- hwglabr::gff_get_attribute(SK1_gff$attributes, 'ID')

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#                           Add all data to package                            #
#                             (as internal data)                               #

# Determine the best compression for the data files
tools::checkRdaFiles('data/')  # Suggests 'bzip2'

setwd('/Users/luis/Google_Drive_NYU/LabShare_Luis/LabWork/Scripts/Rpackages/hwglabr')
devtools::use_data(S288C_conv_midpoint_dist,
                   S288C_div_midpoint_dist,
                   S288C_tand_midpoint_dist,
                   SK1_conv_midpoint_dist,
                   SK1_div_midpoint_dist,
                   SK1_tand_midpoint_dist,
                   S288C_conv,
                   SK1_conv,
                   S288Ccen,
                   SK1cen,
                   SK1rosetta,
                   s288C_gff,
                   SK1_gff,
                   internal = TRUE, overwrite = TRUE, compress = "bzip2")