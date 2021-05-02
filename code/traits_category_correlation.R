traits_category_correlation <- function(data, traits.idx.strat, traits.category, tgt.category){
  dir = "/Users/stephanie/UNC-chapel hill/Spring2021/STOR893/893FinalProject/code"
  source(paste0(dir,"/reorder_cormat.R"))
  source(paste0(dir,"/get_upper_tri.R"))
  
  idx <- which(traits.category == tgt.category )+traits.idx.strat-1 # identify columns in the same category
  tmp <- data[,idx] %>% drop_na() %>% cor() %>% get_upper_tri() %>% melt(na.rm = TRUE)
  ggplot(data = tmp, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile()+
    scale_fill_gradient2(name = "correlation")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
          plot.title = element_text(color="black", size=14, face="bold.italic")) +
    labs(title = paste("Correlation plot for",tgt.category) , x = tgt.category, y = tgt.category) 
}