#https://cran.r-project.org/web/packages/data.tree/vignettes/data.tree.html
library(rpart)
library(data.tree)
# tree framework
tree         <- Node$new("BD-IPMN")
surgery      <- tree$AddChild("surgery") 
surveillance <- tree$AddChild("Surveillance")
comp         <- surgery$AddChild("Complication")
nocomp       <- surgery$AddChild("No complication")
comp_death   <- comp$AddChild("Death")
comp_surv    <- comp$AddChild("Survival")
nocomp_death <- nocomp$AddChild("Death")
nocomp_surv  <- nocomp$AddChild("Survival")
surv_markov  <- surveillance$AddChild("Markov")
nocomp_markov <- nocomp_surv$AddChild("Markov")
comp_markov   <- comp_surv$AddChild("Markov")

print(tree)

popClone <- Clone(tree)
as.data.frame(tree)
ToDataFrameNetwork(tree)



SetGraphStyle(tree, rankdir = "TB")
SetEdgeStyle(tree, arrowhead = "vee", color = "grey35", penwidth = 2)
plot(tree)

