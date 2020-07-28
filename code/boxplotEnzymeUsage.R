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

# Only keep data from relevant GO terms
capUse <- capUse[capUse$protID %in% GO$Entry,]
idx <- match(capUse$protID, GO$Entry)
capUse$GOterm <- gsub('TCA','TCA cycle',GO$system[idx])
colnames(capUse) <- gsub('capUse_','',colnames(capUse))

capUse <- capUse %>% mutate_if(is.numeric, round, digits = 3)
write_delim(capUse,'../../results/enzymeUsage/capUsage.txt',delim = '\t')

capUse <- gather(capUse, 'Condition', 'Usage', 4:8)
capUse$GOterm <- factor(capUse$GOterm, levels=c('Glycolysis','TCA cycle','ETC','Ribosome'))
capUse$Condition <- factor(capUse$Condition, levels=c('CN4','CN22','CN38','CN78','hGR'))

ggplot(capUse, aes(x = Condition, y = Usage, color=GOterm)) +
  geom_boxplot(lwd = 0.35) +
  scale_color_manual(values=c('#CBBBA0','#1D1D1B','#1D71B8','#878787')) +
  facet_grid(. ~ GOterm) +
  labs(x = '', y = 'Capacity usage (%)') + 
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5), text = element_text(size=7), 
        line = element_line(size=0.15), strip.background = element_blank(),
        axis.line = element_line(size=0.15), legend.position='none')
ggsave("selectedGOtermUsage.pdf", width=10, height=4.5, units='cm')
