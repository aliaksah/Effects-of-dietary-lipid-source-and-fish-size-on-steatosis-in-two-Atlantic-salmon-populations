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
tbl_data = read_xlsx(path = "salmon data.xlsx",sheet = 3)[,c(1,3,5:33)]

colnames(tbl_data)



names(tbl_data)[1] = "Size"

targets = names(tbl_data)[3:31]

tbl_data$SID = 2 - as.integer(as.factor(tbl_data$Size)) 
fml = NULL


fml[[1]] = target~ 1 

fml[[2]] = target~ 1 + RapeOil 

fml[[3]] = target~ 1 +  SID 
fml[[4]] = target~ 1  + SID + RapeOil 
fml[[5]] = target~1 + SID  + RapeOil  + I(RapeOil*(SID==0)) + I(RapeOil*(SID==1)) 

tbl_data = tbl_data[1:12,]

reslist = NULL
for (index in 1:length(targets)){
  

  tbl_data$target = tbl_data[[targets[index]]]
  
  
  pcprior = list(prec = list(prior = "pc.prec", param=c(1,0.1)))

  
  
  idx = which(!is.na(tbl_data$target))
  

  
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
  res = cbind(tbl_data$Size[idx],tbl_data$RapeOil[idx], tbl_data$target[idx], summ.inla$linear.predictor[,c(1,3:18)])
  
  names(res)[1:3] = c("Size","RapeOil","target")
  names(res)[-c(1:3)] = paste0("q_",names(res)[-c(1:3)] )
  
  data.residuals = data.frame(resids = unlist(lapply(1:length(idx),FUN = function(i) unlist(lapply(4:20,FUN = function(j) res[i,j]-res[i,3])))))
  
  data.residuals$resids = (data.residuals$resids-mean(data.residuals$resids))/sd(data.residuals$resids)
  data.residuals$preds = unlist(lapply(1:length(idx),FUN = function(i) unlist(lapply(4:20,FUN = function(j) res[i,j]))))
  data.residuals$Size = unlist(lapply(1:length(idx),FUN = function(i) unlist(lapply(4:20,FUN = function(j) res$Size[i]))))
  
  res$rvalue = ks.test(data.residuals$resids, "pnorm")$p.value
  

  plot2 = data.residuals %>% 
    mutate(Fitted = data.residuals$preds,
           Residuals = data.residuals$resids, Size = data.residuals$Size) %>% 
    ggplot(mapping = aes(x = Fitted, y =Residuals, color=Size)) + ggtitle(paste0("residuals for [",index,"]   gene ",targets[index], " and KS p-value ",round(res$rvalue,2))) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0)
  
  
  
  plot3 = data.residuals %>%
    ggplot(aes(x=resids)) + geom_histogram(aes(y = ..density..),
                                           colour = 1, fill = "blue") + ggtitle(paste0("residuals for  [",index,"] accid ",targets[index], " and KS p-value ",round(res$rvalue,2)))+
    geom_density()
  
  
  plot4 = ggplot(data.residuals, aes(sample = resids)) + 
    stat_qq() + stat_qq_line()+ ggtitle(paste0("residuals for [",index,"]  gene ",targets[index], "and KS p-value ",round(res$rvalue,2)))
  
  res$SizeN = paste0("size",res$Size)
  
  plot1 = res %>% ggplot(mapping = aes(x = RapeOil, y = target, color = SizeN, q_0.005quant = q_0.005quant, q_0.995quant = q_0.995quant)) + ggtitle(paste0("               accid ",targets[index]," model ",j," pprob ",round(max(mliks),5))) +  geom_point() +
    geom_line(aes(y = q_mode), size = 1) + expand_limits(x = 0, y = 0) +  geom_ribbon(aes(ymin = q_0.005quant, ymax = q_0.995quant, color = SizeN), alpha = .25) + theme(text = element_text(size = 20))  
  
  grid.arrange(plot1, plot3,plot4, nrow=3)
  
  ggsave(plot = plot1,filename = paste0("plots_digest/digest_",targets[index],"_plot.png"))
  
  reslist[[index]] = list(best = model.inla,probs = mliks, res = res, data.residuals = data.residuals)
  
}


ress = as.data.frame(do.call(rbind, lapply(X = 1:29,function(i){
  probs = sort(reslist[[i]]$probs,decreasing = T)[1:2]
  bf10 = probs[1]/probs[2]
  reslist[[i]]$best$summary.fixed
  
  saveRDS(file = paste0("models_digest/",targets[i],"_digest_best.RDS"),object = (reslist[[i]]$best$summary.fixed))
  
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

View(ress[,c(1:6,19)])

write.csv(ress,"models_digest/final_table.csv")


View(ress)
