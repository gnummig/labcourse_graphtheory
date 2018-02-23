library("ggplot2")
library("gridExtra")
#library("hexbin")

filtered_results <- read.delim("EP_True_filtered_results.txt")
splice_results <- read.delim("EP_True_splice_graph_results.txt")

fr = filtered_results
fr["Flow_IOError"]=abs(fr["In_Flow"]-fr["Out_Flow"])
fr[fr$VertexID<=1,]["Flow_IOError"]=0

orig = fr[fr$GraphKind =='orig',]
comp = fr[fr$GraphKind =='comp',]
res = fr[fr$GraphKind =='res',]

for(i in 1:nrow(res)) {
  row = res[i,]
  if (row["ChimearNode"]>0) {
    orig[orig$GraphID==row$GraphID & orig$VertexID==row$VertexID ,]["ChimearNode"]=1
    comp[comp$GraphID==row$GraphID & comp$VertexID==row$VertexID ,]["ChimearNode"]=1
  }
}

for(i in 1:nrow(comp)) {
  row = comp[i,]
  if (row["ProblemNode"]==1) {
    orig[orig$GraphID==row$GraphID & orig$VertexID==row$VertexID ,]["ProblemNode"]=1
  }
}

for(i in 1:nrow(res)) {
  row = res[i,]
  if (row["ProblemNode"]==1) {
    orig[orig$GraphID==row$GraphID & orig$VertexID==row$VertexID ,]["ProblemNode"]=2
  }
}

sr = splice_results
sr["Flow_IOError"]=abs(sr["In_Flow"]-sr["Out_Flow"])
sr[sr$VertexID<=1,]["Flow_IOError"]=0

orig_sp = sr[sr$GraphKind =='orig',]
comp_sp = sr[sr$GraphKind =='comp',]
res_sp = sr[sr$GraphKind =='res',]

for(i in 1:nrow(res_sp)) {
  row = res_sp[i,]
  if (row["ChimearNode"]>0) {
    orig_sp[orig_sp$GraphID==row$GraphID & orig_sp$VertexID==row$VertexID ,]["ChimearNode"]=1
    comp_sp[comp_sp$GraphID==row$GraphID & comp_sp$VertexID==row$VertexID ,]["ChimearNode"]=1
  }
}

for(i in 1:nrow(comp_sp)) {
  row = comp_sp[i,]
  if (row["ProblemNode"]==1) {
    orig_sp[orig_sp$GraphID==row$GraphID & orig_sp$VertexID==row$VertexID ,]["ProblemNode"]=1
  }
}

for(i in 1:nrow(res)) {
  row = res_sp[i,]
  if (row["ProblemNode"]==1) {
    orig_sp[orig_sp$GraphID==row$GraphID & orig_sp$VertexID==row$VertexID ,]["ProblemNode"]=2
  }
}

#alle Chimäre Knoten:
res_chim = res[sign(res$ChimearNode) ==1,]
res_chim_sp = res_sp[sign(res_sp$ChimearNode) ==1,]
comp_chim = comp[sign(comp$ChimearNode) ==1,]
comp_chim_sp = comp_sp[sign(comp_sp$ChimearNode) ==1,]
#alle ProblemKnoten ohne Chimäre Knoten
res_problem = res[res$ProblemNode ==1 & res$ChimearNode == 0 ,]
res_problem_sp = res_sp[res_sp$ProblemNode ==1 & res_sp$ChimearNode == 0 ,]
comp_problem = comp[comp$ProblemNode ==1 & comp$ChimearNode == 0 ,]
comp_problem_sp = comp_sp[comp_sp$ProblemNode ==1 & comp_sp$ChimearNode == 0 ,]

#Graphen mit beiden Knoten:
#test=merge(chim[c("GraphID", "Vertex_ID", "In_Degree", "Out_Degree", "Centrality")], 
#           problem[c("GraphID","Vertex_ID" ,"In_Degree", "Out_Degree", "Centrality")], 
#           by="GraphID" )

#Statistik
stat_res_chim_filt = sapply(res_chim[c("In_Degree", "Out_Degree","In_Flow", "In_Flow_Std", "Out_Flow", "Out_Flow_Std", "Centrality")], mean)
stat_res_prob_filt = sapply(res_problem[c("In_Degree", "Out_Degree","In_Flow", "In_Flow_Std", "Out_Flow", "Out_Flow_Std", "Centrality")], mean)
stat_res_chim_spl = sapply(res_chim_sp[c("In_Degree", "Out_Degree","In_Flow", "In_Flow_Std", "Out_Flow", "Out_Flow_Std", "Centrality")], mean)
stat_res_prob_spl = sapply(res_problem_sp[c("In_Degree", "Out_Degree","In_Flow", "In_Flow_Std", "Out_Flow", "Out_Flow_Std", "Centrality")], mean)

stat_comp_chim_filt = sapply(comp_chim[c("In_Degree", "Out_Degree","In_Flow", "In_Flow_Std", "Out_Flow", "Out_Flow_Std", "Centrality")], mean)
stat_comp_prob_filt = sapply(comp_problem[c("In_Degree", "Out_Degree","In_Flow", "In_Flow_Std", "Out_Flow", "Out_Flow_Std", "Centrality")], mean)
stat_comp_chim_sp = sapply(comp_chim_sp[c("In_Degree", "Out_Degree","In_Flow", "In_Flow_Std", "Out_Flow", "Out_Flow_Std", "Centrality")], mean)
stat_comp_prob_sp = sapply(comp_problem_sp[c("In_Degree", "Out_Degree","In_Flow", "In_Flow_Std", "Out_Flow", "Out_Flow_Std", "Centrality")], mean)

#Knoten im Originalgraphen die Chimären produzieren
stat_orig_chim_filt = sapply(orig[orig$ChimearNode==1,][c("In_Degree", "Out_Degree","In_Flow", "In_Flow_Std", "Out_Flow", "Out_Flow_Std", "Flow_IOError", "Centrality")], mean)
stat_orig_chim_spl = sapply(orig_sp[orig_sp$ChimearNode==1,][c("In_Degree", "Out_Degree","In_Flow", "In_Flow_Std", "Out_Flow", "Out_Flow_Std", "Flow_IOError", "Centrality")], mean)
#Knoten im Originalgraphen die keine Chimären produzieren
stat_orig_noChim_filt = sapply(orig[orig$ChimearNode==0,][c("In_Degree", "Out_Degree","In_Flow", "In_Flow_Std", "Out_Flow", "Out_Flow_Std", "Flow_IOError", "Centrality")], mean)
stat_orig_noChim_spl = sapply(orig_sp[orig_sp$ChimearNode==0,][c("In_Degree", "Out_Degree","In_Flow", "In_Flow_Std", "Out_Flow", "Out_Flow_Std", "Flow_IOError", "Centrality")], mean)


#### plots ####
#version.labs <- c(`1`="Original Graph", `2`="Compacted Graph", `3`="Resolved Graph", `4`="Version 4.0", `5`="Version 5.0")
my_box_plot <- ggplot(orig, aes(x=factor(ProblemNode), y=In_Degree+Out_Degree, fill=factor(ProblemNode)) ) +
geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(ProblemNode)), alpha=0.9) +
geom_boxplot(alpha = 0.5, show.legend = FALSE) +
#facet_grid(.~, labeller = as_labeller(version.labs)) +
theme(strip.text.x = element_text(size=9, color="black", face="bold"))

print(my_box_plot)

#version.labs <- c(`1`="Original Graph", `2`="Compacted Graph", `3`="Resolved Graph", `4`="Version 4.0", `5`="Version 5.0")
my_box_plot2 <- ggplot(fr, aes(x=factor(GraphKind), y=In_Flow+Out_Flow, fill=factor(GraphKind)) ) +
geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(GraphKind)), alpha=0.9) +
geom_boxplot(alpha = 0.5, show.legend = FALSE) +
#facet_grid(.~, labeller = as_labeller(version.labs)) +
theme(strip.text.x = element_text(size=9, color="black", face="bold"))

print(my_box_plot2)

my_plot <- ggplot() +
  geom_point(data=comp_problem, aes(x=In_Flow_Std,y=Out_Flow_Std)) +
  geom_point(data=comp_chim, aes(x=In_Flow_Std,y=Out_Flow_Std),color='red') +
  #theme(  axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")  ) +
  ylab("CAI / CDS") + xlab("GC content (%) / CDS") +
  geom_point(size=1.5) +
  theme(  axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")  ) 


print(my_plot)

my_plot2 <- ggplot() + 
  geom_point(data=orig[orig$ChimearNode==0,], aes(x=Centrality,y=Flow_IOError)) +
  geom_point(data=orig[orig$ChimearNode==1,], aes(x=Centrality,y=Flow_IOError),color='red') +
  #theme(  axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")  ) +
  ylab("Flow IO/Error") + xlab("Degree Centrality") +
  geom_point(size=1.5) +
  theme(  axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")  ) 


print(my_plot2)
#### scatter plot with density functions
library(gridExtra)
gglegend <- function(x){ 
  tmp <- ggplot_gtable(ggplot_build(x)) 
  leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box") 
  tmp$grobs[[leg]]
}
mydata =  rbind.data.frame(orig[orig$ChimearNode==0,],orig[orig$ChimearNode==1,])
mydata = mydata[which(mydata$Centrality<0.3 & mydata$Flow_IOError< 500),]

#############################
###### density vs node degree 
#############################

#placeholder plot - prints nothing at all
mydummy  <- ggplot(data=mydata,aes(x=Centrality,y=Flow_IOError)) + 
  geom_point(aes(color=as.factor(mydata$ChimearNode))) +
  scale_color_manual(values = c("orange", "purple"),
                     breaks=c(0,1),
                     labels=c("Problem node", "Chimaer node")) +
  #theme(  axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")  ) +
  ylab("Flow IO-Difference") + xlab("Node Degree") +
  #scale_fill_discrete(breaks=c(0,1),
   #                     labels=c("Problem node", "Chimaer node") ) +
  theme(legend.position=c(1,1),
        legend.justification=c(1,1),
        legend.title = element_blank()
        )


empty = gglegend(mydummy)

#scatterplot of x and y variables
scatter <- ggplot(data=mydata,aes(x=Centrality,y=Flow_IOError)) + 
  geom_point(aes(color=as.factor(mydata$ChimearNode))) +
  scale_color_manual(values = c("orange", "purple")) +
  #theme(  axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")  ) +
  ylab("Flow IO-Difference") + xlab("Degree Centrality") +
  theme(legend.position="none")
      #  ,
      #  axis.text=element_text(size=12), 
      #  axis.title=element_text(size=14,face="bold") )
 
#marginal density of x - plot on top
plot_top <- ggplot(mydata, aes(mydata$Flow_IOError, fill=as.factor(mydata$ChimearNode))) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c("orange", "purple")) + 
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(colour="white")
  )



#marginal density of y - plot on the right
plot_right <- ggplot(mydata, aes(mydata$Centrality, fill=as.factor(mydata$ChimearNode))) + 
  geom_density(alpha=.5) + 
  coord_flip() + 
  scale_fill_manual(values = c("orange", "purple")) + 
  theme(legend.position="none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(colour="white")
  )

#arrange the plots together, with appropriate height and width for each row and column
grid.arrange(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

#placeholder plot - prints nothing at all
empty <- ggplot()+geom_point(aes(1,1), colour="white") +
  theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

#scatterplot of x and y variables
my_plot3 <- ggplot() + 
  geom_point(data=orig[orig$ProblemNode==0,], aes(x=Centrality,y=Flow_IOError)) +
  geom_point(data=orig[orig$ProblemNode==2,], aes(x=Centrality,y=Flow_IOError) , color='red') +
  geom_point(data=orig[orig$ProblemNode==1,], aes(x=Centrality,y=Flow_IOError) , color='green') +
  #theme(  axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")  ) +
  ylab("Flow IO/Error") + xlab("Degree Centrality") +
  geom_point(size=1.5) +
  theme(  axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")  ) 

print(my_plot3)
###Ablage###

mean(res[res$GraphID==0,res$In_Degree])
sapply(res[res$GraphID==0,], mean, na.rm=T)

data.frame(res$In_Degree[res$GraphID==0], res$Out_Degree[res$GraphID==0])

#aggregiere über GraphID
aggregate(data.frame(res$GraphID, res$ProblemNode, res$In_Degree, res$Out_Degree), by = list(res$GraphID), FUN=mean, na.rm=TRUE)
#aggregiere über Problem Nodes
aggregate(data.frame(comp$GraphID, comp$ProblemNode, comp$In_Degree, comp$Out_Degree), by = list(comp$ProblemNode), FUN=mean, na.rm=TRUE)

#group-agreggierung
aggregate(.~ProblemNode+GraphID, res[c("GraphID","ProblemNode","In_Degree", "Out_Degree", "Centrality")], mean)
aggregate(.~ProblemNode+GraphKind, fr[c("GraphKind","ProblemNode", "ChimearNode", "In_Degree", "Out_Degree", "Centrality")],
          mean)
aggregate(.~ProblemNode+ChimearNode+GraphKind, fr[c("GraphKind","ProblemNode", "ChimearNode", "In_Degree", "Out_Degree", "Centrality")], mean)
aggregate(In_Degree+Out_Degree~ProblemNode+GraphKind, fr[c("GraphKind","ProblemNode", "ChimearNode", "In_Degree", "Out_Degree", "Centrality")], mean)
