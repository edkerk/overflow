# This script takes the enzyme usage data from ec-models and plots this in various ways
#install.packages("tidyverse") # Install tidyverse if required
library(tidyverse)
# Adjust to the correct directory 
setwd("C:/Work/GitHub/overflow/results/enzymeUsage")

# Load usage information and remove proteins with always zero usage
capUse <- read.delim('enzymeUsages.txt')
capUse <- capUse[,1:8]
capSum <- rowSums(capUse[,4:8])
capUse <- capUse[!capSum==0,]
capUse[,4:8] <- capUse[,4:8]*100

# Plot based on GO term annotation as obtained from Uniprot (first rearrange data)
GO <- read.delim('../../data/selectedAnnotation.txt', stringsAsFactors = F)
glycolysis <- GO$Entry[GO$system == 'Glycolysis']
ribosome <- GO$Entry[GO$system == 'Ribosome']
ETC <- GO$Entry[GO$system == 'ETC']
TCA <- GO$Entry[GO$system == 'TCA']
glycolysis <- capUse[capUse$protID %in% glycolysis,]
glycolysis$subSys = 'Glycolysis'
ribosome <- capUse[capUse$protID %in% ribosome,]
ribosome$subSys = 'Ribosome'
ETC <- capUse[capUse$protID %in% ETC,]
ETC$subSys = 'ETC'
TCA <- capUse[capUse$protID %in% TCA,]
TCA$subSys = 'TCA cycle'

dat <- rbind(glycolysis,ribosome,ETC,TCA)
dat <- unique(dat)
colnames(dat) <- gsub('capUse_','',colnames(dat))
dat <- gather(dat, 'Condition', 'Usage', 4:8)
dat$subSys <- factor(dat$subSys, levels=c('Glycolysis','TCA cycle','ETC','Ribosome'))
dat$Condition <- factor(dat$Condition, levels=c('CN4','CN22','CN38','CN78','hGR'))

ggplot(dat, aes(x = Condition, y = Usage, color=subSys)) +
  geom_boxplot(lwd = 0.35) +
  scale_color_manual(values=c('#CBBBA0','#1D1D1B','#1D71B8','#878787')) +
  facet_grid(. ~ subSys) +
  labs(x = '', y = 'Capacity usage (%)') + 
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5), text = element_text(size=7), 
        line = element_line(size=0.15), strip.background = element_blank(),
        axis.line = element_line(size=0.15), legend.position='none')
ggsave("selectedGOtermUsage.pdf", width=10, height=4.5, units='cm')
