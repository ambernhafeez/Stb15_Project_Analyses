#################
# Drawing genes #
#################

# code is adapted from https://dbsloan.github.io/TS2019/exercises/r_figure_drawing.html
# this code is intended to draw the intron/exon structure of homologs with information extracted from GFF files

# input file should contain information from the GFF such as feature type (exon, gene), start, end, strand, etc for each gene of interest
genes <- read.csv("gene_annotations.txt", sep="\t")
genes$Start <- as.numeric(as.character(genes$Start))
genes$End <- as.numeric(as.character(genes$End))
genes$Feature <- as.factor(genes$Feature)

# calculate length of each annotation
genes$length=genes$End-genes$Start

# function to calculate relative start and end positions
getrelcoords <- function(df, query) {
  df[which(df$HomologQuery==query),"relstart"] = df[which(df$HomologQuery==query), "Start"] - df[which(df$HomologQuery==query & df$Feature=="gene"),"Start"]
  df[which(df$HomologQuery==query),"relend"] = df[which(df$HomologQuery==query), "End"] - df[which(df$HomologQuery==query & df$Feature=="gene"),"Start"]
  return(df)
}

# run getrelcoords function for each unique 'HomologQuery' (annotation)
for (queryseq in unique(genes$HomologQuery)) {
  genes = getrelcoords(df=genes, query=queryseq)
}

# read in a file with phylogenetic information if desired including a position in the clade/tree for each gene
# genes can then be arranged to line up with a tree later
cladeproteins=read.csv("clade_info.csv")
# in this case both dataframes have a gene column which can be used to merge theem
gen=merge(genes,cladeproteins,by="Gene")

head(cladeproteins)

# function to draw polygons
drawannotation <- function(i, df, colour="grey"){
  arrowlen=0.004
  boxsize=0.005
  arrowStart=max(df[i,'prelstart'], df[i,'prelend']- arrowlen)
  polygon(x=c(df[i,'prelstart'],
              arrowStart,
              df[i,'prelend'],
              arrowStart,
              df[i,'prelstart']),
          y=c(df[i,'height'] - boxsize, 
              df[i,'height'] - boxsize, 
              df[i,'height'], 
              df[i,'height'] + boxsize, 
              df[i,'height'] + boxsize),
          col=colour)
}

# separate out the different features into their own dataframes
gen.cds <- gen[which(gen$Feature=="CDS"),]
gen.exon <- gen[which(gen$Feature=="exon"),]
gen.gene <- gen[which(gen$Feature=="gene"),]

#######################################
# plot one line for each row of that contains a gene annotation
#######################################
pdf("gene_structures.pdf", width=10, height=20)
plot.new()

  # plot line of gene length
for(i in 1:nrow(gen.gene)) {
  # draw the line for gene length
  segments(0, gen.gene[i, 'height'], gen.gene[i, 'plength'], gen.gene[i, 'height'], lwd=0.5)
  # if the gene is on the + strand draw a +, else if it's on the - strand draw a -
  if(gen.gene[i,'Strand'] == "+"){
    text(x=gen.gene[i, 'prelstart'] - 0.01, y= gen.gene[i, 'height'], labels= "+")
  } else{
    text(x=gen.gene[i, 'prelstart'] - 0.01, y= gen.gene[i, 'height'], labels= "-")
  }
  # add clade position number
  text(x=gen.gene[i, 'prelstart'] - 0.03, y=gen.gene[i, 'height'], labels=gen.gene[i, 'CladePos'])
}
# draw exons using the drawannotation() function
for(i in 1:nrow(gen.exon)) {
  drawannotation(df=gen.exon, i=i, colour="black")
}
dev.off()