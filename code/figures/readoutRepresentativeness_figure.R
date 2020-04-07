library(tidyverse)
library(CellNOptR)
work_dir = "~/Downloads/"
setwd(work_dir)
node_table = read_csv("NodePropertiesCombiMS.csv")
View(node_table)
midas_annotation = CNOlist("CH003.csv")

readouts = colnames(midas_annotation@signals$`0`)

#add a new column presence/absence of each node in readouts list
node_table = node_table %>% mutate(readout = name %in% readouts)

node_table %>% 
  ggplot(aes(x = readout, y = Indegree)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1)
ggsave("readoutVsOtherNodeIndegree.pdf")

node_table %>% 
  ggplot(aes(x = readout, y = Outdegree)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1)
ggsave("readoutVsOtherNodeOutdegree.pdf")

node_table %>% 
  ggplot(aes(x = Indegree, y = Outdegree, col = readout)) + 
  geom_point()




#### PATHWAY ANALYSIS

# node_table was modified so that names would overlap with gene sets
node_table <- read_tsv('./NodePropertiesCombiMS_mod.csv')
node_table = node_table %>% mutate(readout = name %in% readouts)

# download the list of pathways from MSigDB

library(GSEABase)
genesets_temp = getGmt(con=paste0('~/Downloads/c2.all.v7.0.symbols.gmt'))
genesets_temp = unlist(genesets_temp)
gene_to_term = list()
for (geneset in genesets_temp){
  temp <- geneIds(geneset)
  temp2 <- setName(geneset)
  temp3 <- as.data.frame(cbind(temp, rep(temp2, length(temp))))
  names(temp3) <- c("gene","term")
  gene_to_term[[(length(gene_to_term)+1)]] = temp3
}
gene_to_term <- bind_rows(gene_to_term) %>% as_tibble()

# 1. Interferon Response 
# 2. B-cell receptor
# 3. T-cell receptro 
# 4. Survival (MAPK, etc) 
# 5. Apoptosis 
# 6. Inflammation (X) - hmmmm
# 7. Lipid signaling (X) - okay
# 8. Innate immunity
# 9. Multi-drug response (X) - okay

.f <- function(x, prop='ClosenessCentrality'){
  node_table %>% 
    left_join(gene_to_term %>% 
                filter(term==x) %>% 
                filter(gene%in% node_table$name) %>% 
                transmute(name=gene, term=term)) %>% 
    filter(!is.na(term)) %>% 
    ggplot(aes(x=readout, y=!!sym(prop))) + 
    geom_boxplot()
}

.f <- function(x, prop='ClosenessCentrality'){
  node_table %>% 
    left_join(gene_to_term %>% 
                filter(term==x) %>% 
                filter(gene%in% node_table$name) %>% 
                transmute(name=gene, term=term)) %>% 
    filter(!is.na(term))
}

# these are the pathways we want to have a look at
terms <- c('Interferon'='REACTOME_INTERFERON_SIGNALING',
           "B-cell"='KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY',
           "T-cell"='KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY',
           "Surival"='BIOCARTA_MAPK_PATHWAY',
           "Apoptosis"='KEGG_APOPTOSIS', 
           "Innate immune system"='REACTOME_INNATE_IMMUNE_SYSTEM', 
           'Inflammation'="REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM")
# extract the information for all pathways and save in a list
temp <- map(terms, .f)

# make a single dataframe out of this list
temp <- bind_rows(temp, .id = 'ID')

g <- temp %>% 
  ggplot(aes(x=ID, y=ClosenessCentrality, fill=readout)) + 
  geom_boxplot() + 
    xlab('') + ylab('Closeness Centrality') + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.text.x=element_text(angle = 45, hjust = 1)) + 
    scale_fill_manual(values=c('#009F4D', "#307FE2"), name='', 
                      labels=c('not measured', 'measured'))
ggsave(g, filename = '~/Downloads/review_figure_2.pdf', useDingbats=FALSE,
       width = 8, height = 4)
