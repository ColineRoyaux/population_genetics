#R4.2.2
##Dependencies
library(ape) #5.7-1
library(pegas) #1.3
library(grDevices) #4.2.2

##Parameters entry
#HaploNetInput Table with haplotype, their abundance and the sequence of the haplotype
hapseq_filename <- "" 
#File with matrix of haplotype frequencies for each location
hapfreq_filename <- ""
#Directory where you want to save outputs
dir_out <- ""
#File pattern for haplotype network representation
pat_file_save <- paste0("haplo_net_Bsib", Sys.Date())#, Sys.time())
#Table separator and decimal characters used
sep <- "\t" 
dec <- "."
#Columns hapseq
col_seq <- 5 #sequences
col_hap <- 1 #Haplotype
col_abd <- 2 #Frequence
#Location labels optional
lab_legend <- NULL 
#Haplotype names pattern
pat_haplonames <- "B_sib_COI"

#Epsilon : variable to search for median vectors
E <- 2

##Steps
### Create saving dir
dir_save <- paste0(dir_out, Sys.Date(), "_haplonet_out/")
dir.create(dir_save, recursive = TRUE)
cat(capture.output(ls.str()), sep = "\n", file = paste0(dir_save, pat_file_save, "haplo_net_init_param.txt")) #Get initial parameters

###Load data
hapseq_tab <- read.table(hapseq_filename, sep = sep, dec = dec, header = FALSE)
hapfreq_tab <- read.table(hapfreq_filename, sep = sep, dec = dec, header = TRUE, check.names = FALSE)
###Transform haplotype frequencies table into matrix
hapfreq_tab[is.na(hapfreq_tab)] <- 0 #Replace each NA in frequency matrix by a 0
hapfreq_mat <- as.matrix(hapfreq_tab[, -1]) #Discard first column containing haplotype names and create a matrix
rownames(hapfreq_mat) <- gsub(pat_haplonames, "h", hapfreq_tab[, 1]) #allocate haplotype identifiers to rowames
###Transform table to fit it in DNAbin class
hapseqlow_tab <- t(sapply(strsplit(hapseq_tab[, col_seq],""), tolower)) #Separate one pb by column and transform into lowcase letters to get into a proper format for DNAbin conversion
n_hap <- nrow(hapseqlow_tab) #count haplotype number
rownames(hapseqlow_tab) <- gsub(pat_haplonames, "h", hapseq_tab[, col_hap]) #allocate haplotype identifiers to rowames
#Pegas doesn't handle gaps, Trim the sequences to remove any column containing gaps
tab_hapseq <- as.data.frame(hapseqlow_tab) == "-"
if(any(tab_hapseq)){
  del_pos <- which(apply(tab_hapseq, 2, any))
  cat("\nGaps have been found in the haplotype alignment, as pegas doesn't handle gaps all positions where a gap is found in at least one sequence are removed. \nRemoved following bases for all haplotypes:\n",
      paste0("- ", del_pos, "\n"))
  hapseqlow_tab <- hapseqlow_tab[, -del_pos]
}

hapseq_DNA <- ape::as.DNAbin(hapseqlow_tab) #Convert to DNAbin
#### MJ Network ####

###Create halotype networks
mj_net <- pegas::mjn(hapseq_DNA, prefix = "mv_", epsilon = E)
lab_mj <- attr(mj_net, "labels")
#rmst_mj_net <- pegas::rmst(cophenetic(mj_net), B = 100)
#mj_net <- rmst_mj_net
n_med_vec <- length(grep("mv_", labels(mj_net))) #Count the number of median vectors ("false" haplotype) computed
freq_mat <- rbind(cbind(hapfreq_mat,
                        med_vec = rep(0, n_hap)),
                  cbind(matrix(rep(0, n_med_vec * ncol(hapfreq_mat)),
                               nrow = n_med_vec,
                               ncol = ncol(hapfreq_mat)),
                        med_vec = rep(1, n_med_vec))) #Add median vectors to the frequency matrix => a column dedicated as a "location" to have it represented in a different color
rownames(freq_mat) <- lab_mj
freq_mat_od <- freq_mat[attr(mj_net, "labels"),] #Match order of freq mat and network
#rownames(freq_mat) <- gsub("mv_[0-9]+", "", rownames(freq_mat)) #Remove the names of median vectors
attr(mj_net, "class") <- "haploNet" #Change the class of the network object because else, colors in pie charts doesn't display properly
attr(mj_net, "labels") <- gsub("mv_[0-9]+", "", labels(mj_net)) #Remove the names of median vectors for it to not display
hap_sizes <- apply(freq_mat_od, 1, sum) + 7 #Set sizes for haplotypes
hap_sizes[grep("mv_", names(hap_sizes))] <- 2 #Set sizes for median vectors
pal <- c(grDevices::hcl.colors(ncol(hapfreq_mat), palette = "Set2"), "black")
##Plot OPTIONS
pegas::setHaploNetOptions(pie.colors.function = pal, pie.inner.segments.color = "black", scale.ratio = 6)

if(is.null(lab_legend)){
  lab_legend <- colnames(hapfreq_tab)[-1]
}

##First plot output without alternative paths
svg(paste0(dir_save, pat_file_save, "_mj.svg"), width = 15, height = 7)
plot(mj_net, labels = TRUE, size = hap_sizes, shape = c(rep("circles", n_hap), rep("diamonds", n_med_vec)),
     pie = freq_mat_od, legend = c(20, 16), threshold = 0)#, show.mutation = 1, cex = 1, font = 2, fast = FALSE, col.link = "black", lwd = 1, bg = NULL)
legend(x = 20, y = 10, lab_legend, fill = pal[-length(pal)], cex = 1.25, xpd = TRUE, title = "Haplotypes found in localities:") #colnames(freq_mat), fill = pal, cex = 1.25, xpd = TRUE, title = "Localities")
# rect(xleft = 100, ybottom = 29, xright = 130, ytop = 33, col = "white", border = "white")
#We added 7 to the size of the rounds for better display so we have to correct the legend
# text(x = 103, y = 30, "Haplotype abundance:   0      8    20")
dev.off()

##Second plot output with all alternate paths showing
svg(paste0(dir_save, "alt_", pat_file_save, "_mj.svg"), width = 15, height = 7)
plot(mj_net, labels = TRUE, size = hap_sizes, shape = c(rep("circles", n_hap), rep("diamonds", n_med_vec)),
     pie = freq_mat_od, legend = c(40, 40), threshold = c(1, 1e4))#, show.mutation = 1, cex = 1, font = 2, fast = FALSE, col.link = "black", lwd = 1, bg = NULL)
legend(x = 40, y = 25, lab_legend, fill = pal[-length(pal)], cex = 1.25, xpd = TRUE, title = "Haplotypes found in localities:") #colnames(freq_mat), fill = pal, cex = 1.25, xpd = TRUE, title = "Localities")
# rect(xleft = 100, ybottom = 29, xright = 130, ytop = 33, col = "white", border = "white")
# #We added 7 to the size of the rounds for better display so we have to correct the legend
# text(x = 103, y = 30, "Haplotype abundance:   0      8    20")
dev.off()


#### Parsimonious network ####
haploseq <- pegas::haplotype(hapseq_DNA)
d <- ape::dist.dna(haploseq, "N")

#Haplotype network computation
haplo_net <- pegas::haploNet(haploseq, d = d)
n_med_vec <- 0

freq_mat <- rbind(cbind(hapfreq_mat,
                        med_vec = rep(0, n_hap)),
                  cbind(matrix(rep(0, n_med_vec * ncol(hapfreq_mat)),
                               nrow = n_med_vec,
                               ncol = ncol(hapfreq_mat)),
                        med_vec = rep(1, n_med_vec))) #Add median vectors to the frequency matrix => a column dedicated as a "location" to have it represented in a different color
attr(haplo_net, "labels") <- rownames(freq_mat)

freq_mat_od <- freq_mat[attr(haplo_net, "labels"),] #Match order of freq mat and network
#rownames(freq_mat) <- gsub("mv_[0-9]+", "", rownames(freq_mat)) #Remove the names of median vectors
attr(haplo_net, "class") <- "haploNet" #Change the class of the network object because else, colors in pie charts doesn't display properly
attr(haplo_net, "labels") <- gsub("mv_[0-9]+", "", labels(haplo_net)) #Remove the names of median vectors for it to not display
hap_sizes <- apply(freq_mat_od, 1, sum) + 7 #Set sizes for haplotypes
hap_sizes[grep("mv_", names(hap_sizes))] <- 2 #Set sizes for median vectors
pal <- c(grDevices::hcl.colors(ncol(hapfreq_mat), palette = "Set2"), "black")
##Plot OPTIONS
pegas::setHaploNetOptions(pie.colors.function = pal, pie.inner.segments.color = "black", scale.ratio = 6)

if(is.null(lab_legend)){
  lab_legend <- colnames(hapfreq_tab)[-1]
}

##First plot output without alternative paths
svg(paste0(dir_save, pat_file_save, "_pars.svg"), width = 15, height = 7)
plot(haplo_net, labels = TRUE, size = hap_sizes, shape = c(rep("circles", n_hap), rep("diamonds", n_med_vec)),
     pie = freq_mat_od, legend = c(40, 40), threshold = 0)#, show.mutation = 1, cex = 1, font = 2, fast = FALSE, col.link = "black", lwd = 1, bg = NULL)
legend(x = 40, y = 25, lab_legend, fill = pal[-length(pal)], cex = 1.25, xpd = TRUE, title = "Haplotypes found in localities:") #colnames(freq_mat), fill = pal, cex = 1.25, xpd = TRUE, title = "Localities")
# rect(xleft = 100, ybottom = 29, xright = 130, ytop = 33, col = "white", border = "white")
#We added 7 to the size of the rounds for better display so we have to correct the legend
# text(x = 103, y = 30, "Haplotype abundance:   0      8    20")
dev.off()

##Second plot output with all alternate paths showing
svg(paste0(dir_save, "alt_", pat_file_save, "_pars.svg"), width = 15, height = 7)
plot(haplo_net, labels = TRUE, size = hap_sizes, shape = c(rep("circles", n_hap), rep("diamonds", n_med_vec)),
     pie = freq_mat_od, legend = c(40, 40), threshold = c(1, 1e4))#, show.mutation = 1, cex = 1, font = 2, fast = FALSE, col.link = "black", lwd = 1, bg = NULL)
legend(x = 40, y = 25, lab_legend, fill = pal[-length(pal)], cex = 1.25, xpd = TRUE, title = "Haplotypes found in localities:") #colnames(freq_mat), fill = pal, cex = 1.25, xpd = TRUE, title = "Localities")
# rect(xleft = 100, ybottom = 29, xright = 130, ytop = 33, col = "white", border = "white")
# #We added 7 to the size of the rounds for better display so we have to correct the legend
# text(x = 103, y = 30, "Haplotype abundance:   0      8    20")
dev.off()



#Reset plot options to default

pegas::setHaploNetOptions(labels = TRUE, labels.cex = 1, labels.font = 2, link.color = "black", link.type = 1, link.type.alt = 2, link.width = 1,
                          link.width.alt = 1, haplotype.inner.color = "white", haplotype.outer.color = "black", mutations.cex = 1, mutations.font = 1,
                          mutations.frame.background = "#0000FF4D", mutations.frame.border = "black", mutations.text.color = 1, mutations.arrow.color = "black",
                          mutations.arrow.type = "triangle", mutations.sequence.color = "#BFBFBF4D", mutations.sequence.end = "round",
                          mutations.sequence.length = 0.3,  mutations.sequence.width = 5,  pie.inner.segments.color = "black",
                          pie.colors.function = rainbow, scale.ratio = 1, show.mutation = 1)
dev.off()


##### PARAM BOECKELLA SPI #### 

##First plot output without alternative paths
# svg(paste0(dir_save, pat_file_save, ".svg"), width = 15, height = 7)
# plot(mj_net, labels = TRUE, size = hap_sizes, shape = c(rep("circles", n_hap), rep("diamonds", n_med_vec)),
#      pie = freq_mat_od, legend = c(100, 46), threshold = 0)#, show.mutation = 1, cex = 1, font = 2, fast = FALSE, col.link = "black", lwd = 1, bg = NULL)
# legend(x = 100, y = 28, lab_legend, fill = pal[-length(pal)], cex = 1.25, xpd = TRUE, title = "Haplotypes found in localities:") #colnames(freq_mat), fill = pal, cex = 1.25, xpd = TRUE, title = "Localities")
# rect(xleft = 100, ybottom = 29, xright = 130, ytop = 33, col = "white", border = "white")
# #We added 7 to the size of the rounds for better display so we have to correct the legend
# text(x = 103, y = 30, "Haplotype abundance:   0      8    20")
# dev.off()
# 
# ##Second plot output with all alternate paths showing
# svg(paste0(dir_save, "alt_", pat_file_save, ".svg"), width = 15, height = 7)
# plot(mj_net, labels = TRUE, size = hap_sizes, shape = c(rep("circles", n_hap), rep("diamonds", n_med_vec)),
#      pie = freq_mat_od, legend = c(100, 46), threshold = c(1, 1e4))#, show.mutation = 1, cex = 1, font = 2, fast = FALSE, col.link = "black", lwd = 1, bg = NULL)
# legend(x = 100, y = 28, lab_legend, fill = pal[-length(pal)], cex = 1.25, xpd = TRUE, title = "Haplotypes found in localities:") #colnames(freq_mat), fill = pal, cex = 1.25, xpd = TRUE, title = "Localities")
# rect(xleft = 100, ybottom = 29, xright = 130, ytop = 33, col = "white", border = "white")
# #We added 7 to the size of the rounds for better display so we have to correct the legend
# text(x = 103, y = 30, "Haplotype abundance:   0      8    20")
# 
# 
# dev.off()



