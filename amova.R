#R4.2.2
##Dependencies
library(ade4) #1.7-22

##Parameters entry
haploXpop_filename <- "" #AMOVAinput Haplotype X Population matrix
dist_filename <- "" #AMOVAinput Distance Haplotype matrix
sep <- "\t" #Table separator
dec <- "." #Decimal separator
out_dir <- "" #Output directory

##Steps
### Save init param
cat(capture.output(ls.str()), sep = "\n", file = paste0(out_dir, Sys.Date(),"_Amova_init_param.txt")) #Get initial parameters

###Haplotypes x populations abundance matrix
tab_haploXpop <- read.table(file = haploXpop_filename, sep = sep, dec = dec, header = TRUE)
tab_haploXpop[is.na(tab_haploXpop)] <- 0
tab_haploXpop_c <- tab_haploXpop[,-1]

###Distance Matrix
# Open distance matrix file
tab_dist_haplo <- read.table(file = dist_filename, sep = sep, dec = dec)
dist_haplo <- as.dist(tab_dist_haplo[-1,-1])

# Create Euclidean distance matrix
euc_dist_haplo <- sqrt(dist_haplo)

# Verif if it is euclidean
ade4::is.euclid(euc_dist_haplo)
#TRUE

###Structure table (optional)
#tab_structure <- read.table(file = structure_filename, sep = sep, dec = dec)

###AMOVA TIME
# Perform the amova
Amv <- ade4::amova(tab_haploXpop_c, euc_dist_haplo)
cat(capture.output(Amv), file = paste0(out_dir, Sys.Date(), "_Amova_out.txt"), sep = "\n")

# Test its significance through random test (permutation of matrix)
set.seed(1997)
Amvsignif <- ade4::randtest(Amv, nrepet = 999, alter = "two-sided")

cat("\nSignificance of amova through random permutation: \n", file = paste0(out_dir, Sys.Date(), "_Amova_out.txt"), sep = "\n", append = TRUE)
cat(capture.output(Amvsignif), file = paste0(out_dir, Sys.Date(), "_Amova_out.txt"), sep = "\n", append = TRUE)
png(paste0(out_dir, Sys.Date(), "_sim_Amova_out.png"), width = 600, height = 500)
plot(Amvsignif)
dev.off()

