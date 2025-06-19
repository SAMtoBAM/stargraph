###this is a modifiable Rscript for generating Starship alignments
##the basic skeleton will remain the same and the same features can be swapped in for each variable

options(warn = -1)

##only need gggenomes
suppressMessages(library(IRanges))
suppressMessages(library(gggenomes))
suppressMessages(library(ggnewscale))

##need four features

##First feature
##first a bed file with just the length and size of the full alignment regions
bed=read.csv("PATHTOOUTPUT/2.HGT_candidates/alignments/ELEMENT/ELEMENT.CANDIDATEGENOME2.aligned_regions.bed", sep='\t', header=T)
##just modify some headers for downstream handling
##make another header, copy of contig, but called seq_id
bed$seq_id = bed$contig
##also length using the end position
#bed$length = bed$end
bed$length = bed$end - bed$start
##add a column with label to be used in the plot (showing the actual region being aligned)
bed$label= paste(bed$contig,":",bed$start,"-",bed$end, sep = "")

##second feature
##a bed file of just the Starship-like region coordinates i.e. the above bed file without the flanking regions
SLRbed=read.csv("PATHTOOUTPUT/2.HGT_candidates/alignments/ELEMENT/ELEMENT.element.bed", sep='\t', header=T)
SLRbed$seq_id = SLRbed$contig
SLRbed$length = SLRbed$end-SLRbed$start

##third feature
##a bed file with the genes annotated within the Starship-like regions (coordinates modified due to the flanking regions added)
genes=read.csv("PATHTOOUTPUT/2.HGT_candidates/genes.bed", sep='\t', header=T)
genes$seq_id = genes$contig
genes$length= genes$end-genes$start
genes$strand=genes$sense

##fourth feature
##the nucmer all-v-all alignment converted to paf format
links=suppressMessages(suppressWarnings(read_links("PATHTOOUTPUT/2.HGT_candidates/alignments/ELEMENT/ELEMENT.CANDIDATEGENOME2.contigs.delta.paf")))


##the actual plot
## only selecting alignments greater than 1kb and with greater then 80% identity)
##save plot as variable so can save it
plot=suppressMessages(suppressWarnings(print(gggenomes(genes=genes, seqs=bed, feat=SLRbed, links=subset(links, map_length > 1000 & map_match/map_length > 0.8 & seq_id != seq_id2), adjacent_only = T) %>%
  gggenomes::sync() %>%
  gggenomes::pick() %>%
  gggenomes::flip() +
  geom_link(aes(fill=((map_match/map_length)*100)) ,colour="black", alpha=0.5, offset = 0.05, size=0.1 )+
  scale_fill_gradientn(colours=c("grey100","grey75", "grey50"), name ="Identity (%)", labels=c(80,90,100), breaks=c(80,90,100), limits = c(80, 100))+
  new_scale_fill()+
  geom_seq(linewidth = 0.5)+
  geom_feat(color="red", alpha=.6, linewidth=3)+
  geom_gene(aes(fill=label), stroke=0.1, colour="black", shape = 3)+
  geom_seq_label(aes(label=label))+
  geom_seq_label(aes(label=tag), nudge_y = -.25)+
  geom_gene_tag(aes(label=label), size = 2, nudge_y=0.1, check_overlap = FALSE)+
  scale_fill_manual(values = c("red","blue","lightblue"), breaks=c("tyrR","myb", "duf3723"), name = NULL)+
  theme(legend.position="top", legend.box = "horizontal"))))

##same plot but now allowing for all alignment s against all others
plot2=suppressMessages(suppressWarnings(print(gggenomes(genes=genes, seqs=bed, feat=SLRbed, links=subset(links, map_length > 1000 & map_match/map_length > 0.8 & seq_id != seq_id2), adjacent_only = F) %>%
  gggenomes::sync() %>%
  gggenomes::pick() %>%
  gggenomes::flip() +
  geom_link(aes(fill=((map_match/map_length)*100)) ,colour="black", alpha=0.5, offset = 0.05, size=0.1 )+
  scale_fill_gradientn(colours=c("grey100","grey75", "grey50"), name ="Identity (%)", labels=c(80,90,100), breaks=c(80,90,100), limits = c(80, 100))+
  new_scale_fill()+
  geom_seq(linewidth = 0.5)+
  geom_feat(color="red", alpha=.6, linewidth=3)+
  geom_gene(aes(fill=label), stroke=0.1, colour="black", shape = 3)+
  geom_seq_label(aes(label=label))+
  geom_seq_label(aes(label=tag), nudge_y = -.25)+
  geom_gene_tag(aes(label=label), size = 2, nudge_y=0.1, check_overlap = FALSE)+
  scale_fill_manual(values = c("red","blue","lightblue"), breaks=c("tyrR","myb", "duf3723"), name = NULL)+
  theme(legend.position="top", legend.box = "horizontal"))))

widthFrac = max(bed$length+100000) / 20000
#heightFrac = nrow(regionSeqs) / 2 # for default cases
heightFrac = nrow(bed) / 0.5 # for when OG ids are included as gene_name

suppressMessages(suppressWarnings(ggsave("PATHTOOUTPUT/2.HGT_candidates/alignments/ELEMENT/ELEMENT.CANDIDATEGENOME2.adjacent.png", plot = plot, units = "in", height = heightFrac, width = widthFrac, limitsize = FALSE)))
suppressMessages(suppressWarnings(ggsave("PATHTOOUTPUT/2.HGT_candidates/alignments/ELEMENT/ELEMENT.CANDIDATEGENOME2.adjacent.svg", plot = plot, units = "in", height = heightFrac, width = widthFrac, limitsize = FALSE)))

suppressMessages(suppressWarnings(ggsave("PATHTOOUTPUT/2.HGT_candidates/alignments/ELEMENT/ELEMENT.CANDIDATEGENOME2.allvall.png", plot = plot2, units = "in", height = heightFrac, width = widthFrac, limitsize = FALSE)))
suppressMessages(suppressWarnings(ggsave("PATHTOOUTPUT/2.HGT_candidates/alignments/ELEMENT/ELEMENT.CANDIDATEGENOME2.allvall.svg", plot = plot2, units = "in", height = heightFrac, width = widthFrac, limitsize = FALSE)))

