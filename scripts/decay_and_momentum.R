library("ggplot2")
library("tidyverse")
library("ggrepel")
library("patchwork")


options(stringsAsFactors=F)
options(scipen=99999)

mytheme = theme_bw()+ theme(text = element_text(size=12))

catcolors=c("GAM.NPMI"="#c83737",
            "meanval.logval1"="#0063FF",
            "meanval.logval2"="#525252",
            "meanval.logSPRITE"="#ff9700"
)


chrOrder<-paste0("chr",c(1:22,"X","Y"))

#round_any = function(x, accuracy, f=floor){f(x/ accuracy) * accuracy}
#
#robust.scaling = function(x){
#  q=quantile(x, c(0.01, 1), na.rm = T)
#  r=(x - min(q)) / diff(q)
#  return(r)
#}

robust.scaling = function(x, x2, p.min=0.05, p.max=.95){
  #scale values based on quantiles so we max the overlap and exclude outlayers at the tails
  q=quantile(x2, c(p.min, p.max), na.rm = T)
  r=(x - min(q, na.rm = T)) / diff(range(q, na.rm = T))
  return(r)
}




###### get data

read_GAM = function(fn){
  mat= read_tsv(fn)
  lng = mat %>% pivot_longer(-1) %>% tidyr::separate(1, into=c("chrom_x", "start_x", "end_x")) %>%
    tidyr::separate(name, into=c("chrom_y", "start_y", "end_y")) %>%
    dplyr::filter(!is.na(value)) %>%
    mutate(across(c("start_x", "end_x","start_y", "end_y"), as.numeric)) %>%
    mutate(bindist=start_y-start_x) %>% 
    filter(bindist>1)
  return(lng)
}

file_list_GAM = dir("/data/pombo/Sasha/gamtools/H1_hESC_matrices/hg38/Curated/20PMR_60OW/npmi/per_chromosome/at50Kb/", full.names = T)
#contacts.GAM = file_list_GAM %>% map_dfr(read_GAM, .id = "ds")
#write_rds(contacts.GAM, "/data/pombo/christoph/2201_H1_distance_decay/data/contacts.GAM2.rds")

contacts.GAM = read_rds("/data/pombo/christoph/2201_H1_distance_decay/data/contacts.GAM2.rds")

contacts.GAM.filtered.means =
  contacts.GAM %>% filter(bindist>0) %>%
  group_by(chrom_x, bindist) %>% 
  summarize(name="GAM.NPMI", value = mean(value, na.rm = T)) %>%
  arrange(chrom_x, bindist) 


#### SPRITE
read_sprite = function(fn){
  mat= read_tsv(fn, col_names = c("chrom_x", "start_x", "end_x","chrom_y", "start_y", "end_y", "SPRITE"))
  lng = mat %>% 
    #dplyr::filter(!is.na(value)) %>%
    mutate(bindist=start_y-start_x) %>% 
    filter(bindist>0)
  return(lng)
}

file_list_sprite = dir("/data/pombo/Teresa/fanc/SPRITE_4DNFITX6WCRT/", pattern = "*50kb*", full.names = T)
file_list_sprite

contacts.SPRITE = file_list_sprite %>% map_dfr(read_sprite, .id = "ds")


contacts.SPRITE.filtered.means =
  contacts.SPRITE %>% filter(bindist>0) %>%
  group_by(chrom_x, bindist) %>% 
  summarize(meanval.logSPRITE = log(mean(SPRITE, na.rm = T))
  ) %>%
  arrange(chrom_x, bindist) 


#download 4DNFI82R42AD
#
#wget https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/cc44f5a1-9d53-480a-99fb-b9a8573f2f26/4DNFI82R42AD.mcool
#file_list_HiC = dir("/data/pombo/christoph/datasets/4DNFI82R42AD/", full.names = T)


read_HiC = function(fn){
  mat= read_tsv(fn, col_names = c("chrom_x", "start_x", "end_x","chrom_y", "start_y", "end_y", "val1", "val2"))
  lng = mat %>% 
    #dplyr::filter(!is.na(value)) %>%
    mutate(bindist=start_y-start_x) %>% 
    filter(bindist>0)
  return(lng)
}


file_list_HiC = dir("/data/pombo/Teresa/fanc/4DNFI82R42AD/long/", full.names = T)
contacts.HiC = file_list_HiC %>% map_dfr(read_HiC, .id = "ds")
#write_rds(contacts.GAM, "/data/pombo/christoph/2201_H1_distance_decay/data/contacts.GAM2.rds")

hist(log(contacts.HiC$val1))
hist(log(contacts.HiC$val2))


contacts.HiC.filtered.means =
  contacts.HiC %>% filter(bindist>0) %>%
  group_by(chrom_x, bindist) %>% 
  summarize(meanval.logHiC.raw = log(mean(val1, na.rm = T)),
            meanval.logHiC.ICE = log(mean(val2, na.rm = T))
  ) %>%
  arrange(chrom_x, bindist) 
#%>%
#  mutate(slope.val1=meanval.val1-lag(meanval.val1),
#         slope.val2=meanval.val2-lag(meanval.val2))


contacts.filtered.means =
  contacts.HiC.filtered.means %>% pivot_longer(-c(chrom_x, bindist) ) %>% # reshape Hi-C fields
  bind_rows(contacts.GAM.filtered.means) %>%
  bind_rows(contacts.SPRITE.filtered.means %>% pivot_longer(-c(chrom_x, bindist)))


###### compute meand and align data 
PMIN=0
PMAX=1

scaling.point.lower=2*50000
scaling.point.upper=100*50000



contacts.filtered.means = contacts.filtered.means %>% 
  #flag last 100 bins from each chromosome to remove noise at very long range distances
  group_by(name, chrom_x) %>%
  mutate(minpos = rank(bindist),
         maxpos = rank(desc(bindist)),
         value.flagged=ifelse(bindist>=scaling.point.lower & bindist<=scaling.point.upper, value, NA)) %>%
  mutate(value.scaled = robust.scaling(value, value.flagged, PMIN, PMAX)) # scale values [0:1]


decay.plot = 
  ggplot(contacts.filtered.means, aes(x=log10(bindist), y=value.scaled, color=name))+
  geom_line(alpha=0.7)+
  geom_vline(xintercept = log10(scaling.point.lower), linetype=2, alpha=.3)+
  geom_vline(xintercept = log10(scaling.point.upper), linetype=2, alpha=.3)+
  ylab("Relative score \n [min:max]")+
  facet_wrap(.~factor(chrom_x, levels=chrOrder), nrow = 2)
decay.plot



BREAKPOINT.DIST=0.1

contacts.filtered.binned.means = 
  contacts.filtered.means %>% 
  group_by(name, chrom_x) %>%
  mutate(bindist.log10 = log10(bindist), bindist.loggroup = cut_width(bindist.log10, BREAKPOINT.DIST)) %>% # groups of constant width
  group_by(name, chrom_x, bindist.loggroup) %>% filter(bindist.log10==min(bindist.log10)) %>%
  group_by(name, chrom_x) %>% arrange(chrom_x, bindist.log10) %>%
  mutate(absdiff=value.scaled-lag(value.scaled), bindist.grouped=(bindist-lag(bindist)), bindiff.grouped.log10=log10(bindist.grouped), slope=absdiff/bindiff.grouped.log10 )


ksf = function(x,y){
  fit = ksmooth(x, y, 'normal',bandwidth=.3) #
  return(list("x"=fit$x, "y"=fit$y))
}


contacts.filtered.binned.means.smoothed = contacts.filtered.binned.means %>% 
  filter(!is.na(slope)) %>% ungroup() %>%
  group_by(chrom_x, name) %>% arrange(chrom_x, bindist.log10) %>% summarize(x=ksf(bindist.log10, slope)$x, y=ksf(bindist.log10, slope)$y)


# show breakpoints
decay.plot2 = 
  ggplot(contacts.filtered.means, aes(x=bindist, y=value.scaled, color=name))+
  geom_line()+
  geom_point(data=contacts.filtered.binned.means)+
  #scale_x_log10()+
  ylab("Relative score \n [5perc:max]")+
  #facet_wrap(.~factor(chrom_x, levels=paste0("chr",c(1:5) )), nrow = 1, scales="free_y")
  facet_wrap(.~factor(chrom_x, levels=chrOrder), nrow = 2)



momentum.plot=ggplot(contacts.filtered.binned.means, aes(x=bindist.log10, y=slope, color=name))+
  geom_point(alpha=0.4, size=1)+
  #geom_line(alpha=0.4)+
  geom_line(data=contacts.filtered.binned.means.smoothed, aes(x=x, y=y, color=name), size=1)+
  ylab("Slope")+
  #facet_wrap(.~factor(chrom_x, levels=paste0("chr",c(1:5) )), nrow = 1, scales="free_y")+
  facet_wrap(.~factor(chrom_x, levels=chrOrder), nrow = 2)+
  ylim(c(-0.05, 0.015))
momentum.plot


(decay.plot & ggtitle("Distance decay")) / (momentum.plot & ggtitle("Momentum")) & 
  scale_color_manual(labels = c("GAM NPMI", "Hi-C counts (log)", "Hi-C ICE (log)", "SPRITE (log)"),
                     values = catcolors) & xlim(c(4.6,8.4)) & 
  theme_minimal() &
  xlab("Genomic distance (log10)")

# annotation
decay.plot.chrom1 = 
  ggplot(contacts.filtered.means %>% filter(chrom_x=="chr1"), aes(x=log10(bindist), y=value.scaled, color=name))+
  geom_line(alpha=0.7)+
  ylab("Relative score \n [5perc:max]")+
  facet_wrap(.~factor(chrom_x, levels=chrOrder), nrow = 2)
#decay.plot.chrom1

momentum.plot.chrom1=ggplot(contacts.filtered.binned.means %>% filter(chrom_x=="chr1"), aes(x=bindist.log10, y=slope, color=name))+
  geom_point(alpha=0.4, size=1)+
  #geom_line(alpha=0.4)+
  geom_line(data=contacts.filtered.binned.means.smoothed %>% filter(chrom_x=="chr1"), aes(x=x, y=y, color=name), size=1)+
  ylab("Slope")+
  #facet_wrap(.~factor(chrom_x, levels=paste0("chr",c(1:5) )), nrow = 1, scales="free_y")+
  facet_wrap(.~factor(chrom_x, levels=chrOrder), nrow = 2)+
  ylim(c(-0.05, 0.015))
#momentum.plot.chrom1


(decay.plot.chrom1 & ggtitle("Distance decay")) / (momentum.plot.chrom1 & ggtitle("Momentum")) & 
  scale_color_manual(labels = c("GAM NPMI", "Hi-C counts (log)", "Hi-C ICE (log)", "SPRITE (log)"),
                     values = catcolors) & xlim(c(4.6,8.4)) & 
  theme_minimal() &
  xlab("Genomic distance (log10)") &
  
  geom_vline(xintercept = 5, linetype=2, alpha=.3) &
  #geom_vline(xintercept = 6, linetype=2, alpha=.3) &
  geom_vline(xintercept = 7, linetype=2, alpha=.3) &
  geom_vline(xintercept = 7.4, linetype=2, alpha=.3) 
#geom_vline(xintercept = , linetype=2, alpha=.3)



###### check variance in adjacent regions
WIDTH=10*50000

contacts.filtered.binned.selected = contacts.filtered.means %>% ungroup() %>% mutate(distancegroup = 
                                                                                       case_when(#10 bins width
                                                                                         between(bindist, 10^5, 10^5+WIDTH) ~ "5",
                                                                                         between(bindist, 10^6, 10^6+WIDTH) ~ "6",
                                                                                         between(bindist, 10^7, 10^7+WIDTH) ~ "7",
                                                                                         between(bindist, 10^7.4, 10^7.4+WIDTH) ~ "7.4",
                                                                                         between(bindist, 10^8, 10^8+WIDTH) ~ "8"
                                                                                       )
)


contacts.filtered.binned.selected2 = contacts.filtered.binned.selected  %>% filter(!is.na(distancegroup)) %>% 
  group_by(chrom_x, name, distancegroup) %>% 
  summarize(centered = value.scaled-mean(value.scaled)      )

value_spread.plot=ggplot(contacts.filtered.binned.selected2, aes(x=as.factor(distancegroup), y=centered, color=name))+
  geom_boxplot()+
  ylab("Spread (centered mean values)")+
  facet_wrap(factor(chrom_x, levels=c("chr1"))~., ncol =4, scales = "free_y")

value_spread.plot & scale_color_manual(labels = c("GAM NPMI", "Hi-C counts (log)", "Hi-C ICE (log)", "SPRITE (log)"),
                                       values = catcolors) & 
  ggtitle("Spread of mean intensities within 10 adjacent bins") &
  theme_minimal() &
  xlab("Genomic distance (log10)") 



