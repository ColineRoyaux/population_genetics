#R4.2.2
##Dependencies
library(adegenet) #2.1.10
library(mmod) #1.3.3

##Parameters entry
haploindiv_filename <- "" #Diffinput Table with individuals, their location, haplotype and population
sep <- "\t" #Table separator
dec <- "." #Decimal separator
out_dir <- "" #Output directory

#Column containing haplotypes
haplo_col <- 1

#Column containing individuals
indivnames_col <- 2

#Column containing populations
pop_col <- 36

##Steps
### Save init param
cat(capture.output(ls.str()), sep = "\n", file = paste0(out_dir, Sys.Date(),"_difftest_init_param.txt")) #Get initial parameters

##Load data
haploindiv_tab <- read.table(haploindiv_filename, sep = sep, dec = dec, header = FALSE)
hap_tab <- as.data.frame(haploindiv_tab[, haplo_col])
row.names(hap_tab) <- haploindiv_tab[, indivnames_col]
colnames(hap_tab) <- "LocusA"

#Convert dataset
hap_gen <- adegenet::df2genind(hap_tab, sep = sep)

#Create population information
adegenet::pop(hap_gen) <- haploindiv_tab[, pop_col]

###Perform test
dif_t <- mmod::diff_test(hap_gen)
cat(capture.output(dif_t), file = paste0(out_dir, Sys.Date(), "_diff_test_out.txt"), sep = "\n")

