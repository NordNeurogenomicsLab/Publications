suppressPackageStartupMessages(library(ggplot2));suppressPackageStartupMessages(library(ggrepel))
w<-Sys.getenv('JIA_WORK_DIR', unset='/mnt/disks/data1/Jia_GEO_submission/codex_review')
set.seed(17)
for(n in c('VZ','MZ')){
 d<-read.csv(file.path(w,'inputs/DE_tables/From_previous_integrations',paste0(n,'_markers.csv')))
 d$plot_log2FC<- -d$avg_log2FC
 d$direction<-ifelse(d$plot_log2FC>0,'High at E13.5','High at E11.5')
 labels<-head(d[order(d$p_val),],40)
 p<-ggplot(d,aes(plot_log2FC,-log10(pmax(p_val,.Machine$double.xmin)),color=direction))+geom_point(size=1,alpha=.5)+geom_text_repel(data=labels,aes(label=gene_name),size=3,max.overlaps=Inf,seed=17)+scale_color_manual(values=c('High at E11.5'='#F8766D','High at E13.5'='#7F8CCB'))+theme_light()+labs(title=paste('E11.5 vs E13.5',n,'markers'),x='Log2 fold change (E13.5 / E11.5)',y='-Log10 P',color=NULL)+theme(legend.position='bottom')
 ggsave(file.path(w,'figures',paste0(if(n=='VZ') 'Fig3B' else 'Fig4A','_volcano.png')),p,width=7,height=7,dpi=220,bg='white')
 write.csv(d,file.path(w,'audit',paste0(n,'_volcano_data.csv')),row.names=FALSE)
}
