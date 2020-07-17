# This script takes the enzyme usage data from ec-models and plots this in various ways
#install.packages("tidyverse") # Install tidyverse if required
library(tidyverse)
# Adjust to the correct directory 
setwd("C:/Work/GitHub/overflow/results/enzymeUsage")

# Load usage information and remove proteins with always zero usage
capUse <- read.delim('enzymeCapUsages.txt')
capSum <- rowSums(capUse[,3:7])
capUse <- capUse[!capSum==0,]
capUse[,3:7] <- capUse[,3:7]*100

# Plot based on GO term annotation as obtained from Uniprot (first rearrange data)
GO <- read.delim('../../data/selectedAnnotation.txt', stringsAsFactors = F)
glycolysis <- GO$Entry[GO$system == 'Glycolysis']
ribosome <- GO$Entry[GO$system == 'Ribosome']
ETC <- GO$Entry[GO$system == 'ETC']
glycolysis <- capUse[capUse$protID %in% glycolysis,]
glycolysis$subSys = 'Glycolysis'
ribosome <- capUse[capUse$protID %in% ribosome,]
ribosome$subSys = 'Ribosome'
ETC <- capUse[capUse$protID %in% ETC,]
ETC$subSys = 'ETC'

dat <- rbind(glycolysis,ribosome,ETC)
dat <- unique(dat)
dat <- gather(dat, 'Condition', 'Usage', 3:7)
dat$subSys <- factor(dat$subSys, levels=c('Glycolysis','ETC','Ribosome'))

ggplot(dat, aes(x = Condition, y = Usage, color=subSys)) +
  geom_boxplot(lwd = 0.35) +
  scale_color_manual(values=c('#C4C6C6','#64A5A3','#72788D')) +
  facet_grid(. ~ subSys) +
  labs(x = '', y = 'Capacity usage (%)', color='') + 
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5), text = element_text(size=7), 
        line = element_line(size=0.15), strip.background = element_rect(size=0.15),
        axis.line = element_line(size=0.15))
ggsave("selectedGOtermUsage.pdf", width=10, height=4.5, units='cm')
  