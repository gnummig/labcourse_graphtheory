filtered_results <- read.delim("~/labcourse_graphtheory/filtered_results.txt")
splice_results <- read.delim("~/labcourse_graphtheory/splice_results.txt")

fr = filtered_results
fr["Flow_IOError"]=abs(fr["In_Flow"]-fr["Out_Flow"])
fr[fr$VertexID<=1,]["Flow_IOError"]=0

orig = fr[fr$GraphKind =='orig',]
comp = fr[fr$GraphKind =='comp',]
res = fr[fr$GraphKind =='res',]

for(i in 1:nrow(res)) {
  row = res[i,]
  if (row["ChimearNode"]==1) {
    orig[orig$GraphID==row$GraphID & orig$VertexID==row$VertexID ,]["ChimearNode"]=1
    comp[comp$GraphID==row$GraphID & comp$VertexID==row$VertexID ,]["ChimearNode"]=1
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
  if (row["ChimearNode"]==1) {
    orig_sp[orig_sp$GraphID==row$GraphID & orig_sp$VertexID==row$VertexID ,]["ChimearNode"]=1
    comp_sp[comp_sp$GraphID==row$GraphID & comp_sp$VertexID==row$VertexID ,]["ChimearNode"]=1
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
