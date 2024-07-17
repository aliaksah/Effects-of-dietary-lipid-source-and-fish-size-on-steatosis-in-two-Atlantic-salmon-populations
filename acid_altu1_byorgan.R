############################################
### Load packages
###
library(readxl)
library(tidyverse)
library(INLA)
library(gridExtra)   
############################################
### Import data from excel
###
tbl_data <- read_xlsx(path = "Fish and fatty acid results for database Final.xlsx",
                      sheet = 4)[1:48,]

colnames(tbl_data)
colnames(tbl_data)


targets = colnames(tbl_data)[4:19]


tbl_data$OID = as.integer(as.factor(tbl_data$Organ)) - 1
tbl_data$OID[tbl_data$Organ=="Feed"] = 1
tbl_data$OID[tbl_data$Organ=="PI"] = 2
tbl_data$OID[tbl_data$Organ=="LI"] = 3
tbl_data$OID[tbl_data$Organ=="Mes"] = 4
tbl_data$SID = 2 - as.integer(as.factor(tbl_data$Fish)) 
tbl_data$LID = as.integer(tbl_data$Organ=="Liver")
tbl_data$LOID = (tbl_data$LID)*(tbl_data$SID+1)

fml = NULL




fml[[1]] = target~ 1 

fml[[2]] = target~ 1 + RapeOil 

fml[[3]] = target~ 1 + SID

fml[[4]] = target~ 1  + SID + RapeOil 

fml[[5]] = target~ 1   + SID  + RapeOil + I(RapeOil*(SID==0)) + I(RapeOil*(SID==1))  



reslist = NULL
for (index in 1:length(targets)){
  tbl_data$target = tbl_data[[targets[index]]]
  
  resaccid = NULL
  
  for(orgid in 1:4)
  {
    
    org = unique(tbl_data$Organ)[orgid]
    
    pcprior = list(prec = list(prior = "pc.prec", param=c(1,0.1)))
    
    
    idx = 1:dim(tbl_data)[1]
    
    
    
    idx = which(tbl_data$Organ==org)
    
    
    
    model.inla = inla(fml[[1]],family = "gaussian",quantiles=c(0.025, 0.5, 0.975,0.005,0.995,0.05,0.95,0.8,0.2,0.9,0.1, 0.7,0.3,0.6,0.4),data = tbl_data[idx,],control.compute=list(config = TRUE,dic=TRUE, waic=TRUE),
                      control.predictor=list(compute=T,quantiles=c(0.025, 0.5, 0.975,0.005,0.995,0.05,0.95,0.8,0.2,0.9,0.1, 0.7,0.3,0.6,0.4)))
    i = 1
    j = 1
    mliks = model.inla$mlik[[1]]
    for(i in 2:5)
    {
      print(i)
      model.inla.cur = inla(fml[[i]],family = "gaussian",quantiles=c(0.025, 0.5, 0.975,0.005,0.995,0.05,0.95,0.8,0.2,0.9,0.1, 0.7,0.3,0.6,0.4),data = tbl_data[idx,],control.compute=list(config = TRUE,dic=TRUE, waic=TRUE),
                            control.predictor=list(compute=T,quantiles=c(0.025, 0.5, 0.975,0.005,0.995,0.05,0.95,0.8,0.2,0.9,0.1, 0.7,0.3,0.6,0.4)))
      mliks = c(mliks,model.inla.cur$mlik[[1]])
      if(model.inla$mlik[[1]]<model.inla.cur$mlik[[1]])
      {
        j = i
        model.inla =  model.inla.cur
      }
    }
    
    mliks = exp(mliks - min(mliks))/sum(exp(mliks - min(mliks)))
    
    print(round(mliks,2))
    
    
    summ.inla = summary(model.inla)
    
    res = cbind(tbl_data$Organ[idx],tbl_data$RapeOil[idx], tbl_data$target[idx], summ.inla$linear.predictor[,c(1,3:18)],tbl_data$Fish[idx])
    names(res)[c(1:3,length(names(res)))] = c("Organ","RapeOil","Acid_conc","Fish")
    names(res)[4:20] = paste0("q_",names(res)[4:20])
    res$order = j
    
    
    
    
    
    data.residuals = data.frame(resids = unlist(lapply(1:length(idx),FUN = function(i) unlist(lapply(4:20,FUN = function(j) res[i,j]-res$Acid_conc[i])))))
    
    data.residuals$resids = (data.residuals$resids-mean(data.residuals$resids))/sd(data.residuals$resids)
    data.residuals$preds = unlist(lapply(1:length(idx),FUN = function(i) unlist(lapply(4:20,FUN = function(j) res[i,j]))))
    data.residuals$organ = unlist(lapply(1:length(idx),FUN = function(i) unlist(lapply(4:20,FUN = function(j) res$Organ[i]))))
    
    res$rvalue = ks.test(data.residuals$resids, "pnorm")$p.value
    
    resaccid = rbind(resaccid,res)
    
    
    
    plot2 = data.residuals %>% 
      mutate(Fitted = data.residuals$preds,
             Residuals = data.residuals$resids, Organ = data.residuals$organ) %>% 
      ggplot(mapping = aes(x = Fitted, y =Residuals, color=Organ)) + ggtitle(paste0("residuals for [",index,"]   accid ",targets[index], " and KS p-value ",round(res$rvalue,2))) +
      geom_point(size = 3) +
      geom_hline(yintercept = 0)
    
    
    
    plot3 = data.residuals %>%
      ggplot(aes(x=resids)) + geom_histogram(aes(y = ..density..),
                                             colour = 1, fill = "blue") + ggtitle(paste0("residuals for  [",index,"] accid ",targets[index], " and KS p-value ",round(res$rvalue,2)))+
      geom_density()
    
    
    plot4 = ggplot(data.residuals, aes(sample = resids)) + 
      stat_qq() + stat_qq_line()+ ggtitle(paste0("residuals for [",index,"]  accid ",targets[index], "and KS p-value ",round(res$rvalue,2)))
    
    
    grid.arrange(plot3,plot4, nrow=2)
    
    
    reslist[[(index-1)*4+orgid]] = list(best = model.inla,probs = mliks, res = res, data.residuals = data.residuals)
  }
  
  
  resaccid$Legend = paste0(resaccid$Organ,ifelse(resaccid$order>2,paste0(" : ",resaccid$Fish),"")," ",": m",resaccid$order)
  plot1 = resaccid %>% ggplot(mapping = aes(x = RapeOil, y = Acid_conc, color = Legend)) + ggtitle(paste0(targets[index])) +  geom_point() +
    geom_ribbon(aes(ymin = q_0.025quant, ymax = q_0.975quant, fill = Legend, color = NULL), alpha = .25) + expand_limits(x = 0, y = 0)+ #ylim(0,ifelse(targets[index]=="C18_1n9c",390,ifelse(targets[index] %in% c("C16_0","C18_2n6c"),150,80)))+
    geom_line(aes(y = q_mode), size = 1) + theme(text = element_text(size = 20))  
  grid.arrange(plot1, nrow=1)
  
  ggsave(plot = plot1,filename = paste0("plots_acid_altu1/acid_",targets[index],"_plot.png"))
  
  
}

targets = colnames(tbl_data)[4:19]

targets =  unlist(lapply(targets, function(tar) paste0(tar," ",unique(tbl_data$Organ)[1:4])))


ress = as.data.frame(do.call(rbind, lapply(X = 1:length(reslist),function(i){
  probs = sort(reslist[[i]]$probs,decreasing = T)[1:2]
  bf10 = probs[1]/probs[2]
  reslist[[i]]$best$summary.fixed
  
  #saveRDS(file = paste0("models_acid_altu1/",targets[i],"_acid_best.RDS"),object = (reslist[[i]]$best$summary.fixed))
  
  id.best =  which(reslist[[i]]$probs == probs[1])
  
  effects = array(0,15)
  effects[1:3] = unlist(reslist[[i]]$best$summary.fixed[1,c(18,3,5)])
  if(id.best == 2)
    effects[7:9] = unlist(reslist[[i]]$best$summary.fixed[2,c(18,3,5)])
  if(id.best == 4)
    effects[7:9] = unlist(reslist[[i]]$best$summary.fixed[3,c(18,3,5)])
  if(id.best == 3 || id.best  == 4|| id.best  == 5)
    effects[4:6] = unlist(reslist[[i]]$best$summary.fixed[2,c(18,3,5)])
  if(id.best  == 5)
  {  
    effects[10:12] = unlist(reslist[[i]]$best$summary.fixed[4,c(18,3,5)])
    effects[13:15] = unlist(reslist[[i]]$best$summary.fixed[5,c(18,3,5)])
  }
  c(paste0("m",id.best),round(c(probs[1],bf10,effects,reslist[[i]]$res$rvalue[1]),3))
})))



rownames(ress) = targets

colnames(ress) = c("Best","Prob","BF12","B0","B0L","B0U","BS","BSL","BSU","BR","BRL","BRU","BRS0","BRS0L","BRS0U","BRS1","BRS1L","BRS1U","KS")

View(ress[,c(1:9,19)])

View(ress)

write.csv(ress,"models_acid_altu1/final_table.csv")
