rm(list=ls())
library(tidyverse)
library(plotly)
library(reticulate)
library(ade4)
library(reticulate)
library(htmltools)
library('Ternary')
dirName <- dirname(rstudioapi::getSourceEditorContext()$path )
setwd(dirName)#setwd(dir = "/scratch/anissa.el/ImmuneStates/")


MarkersCellsTMENs <- readRDS("./MarkerSCellsTMENs.rds")
### CONSTANT VARIABLES
NBCELLTYPES = 17
CELLTYPES <- c('CD8-T', 'Other immune', 'DC / Mono', 'CD3-T', 'B', 'NK', 'Keratin-positive tumor', 'Tumor','CD4-T', 'Mesenchymal-like', 'Macrophages', 'Endothelial', 'Tregs', 'Unidentified', 'DC', 'Mono / Neu','Neutrophils')
NB_TMENs <- 4 
eps <- 10^-4 ## small variation we add to alpha values in case they are null to avoid pbs in the fitting

### FUNCTION TO OPTIMIZE
func <-function(data, pars){

  with(data, sum((value - (pars[NB_TMENs+1]* apply(cbind(arch1,arch2,arch3,arch4),1,function (x) prod(diag(outer(x,pars[1:NB_TMENs],FUN="^"))))))^2))#with(data, sum((value - (pars[NB_TMENs+1]* prod(diag(outer(c(arch1,arch2,arch3,arch4),pars[1:NB_TMENs],FUN="^")))))^2)) #(as.vector(x_cm[,"value"]) - x_cap)^2
}
print("check 1")
### FUNCTION TO GET INPUT DATA: INTENSITY FOR MARKER M FROM CELL TYPE C + TMENs position
get_tmens_coord_marker <- function(MarkersCellsTMENs,cellType, marker,eps){
  MarkersCellsTMENS.ct <- MarkersCellsTMENs%>%
    filter(cell_type == cellType & marker == marker)#%>%
  #pivot_wider(id_cols = c("site_id","patient_id","arch1","arch2","arch3","arch4"),
  #            names_from = "marker",
  #            values_from = "value")
  # machine precision is 2.220446e-16: least positive number x such as 1+ x!  = 1
  # 0 is < mchine precision
  # each time we find a nb that is null, assign this value to avoid Nan, etc...
  MarkersCellsTMENS.ct <- MarkersCellsTMENS.ct  %>%mutate(arch1=ifelse(arch1<.Machine$double.eps,eps,arch1))%>%
    mutate(arch2=ifelse(arch2 <.Machine$double.eps,eps,arch2))%>%
    mutate(arch3=ifelse(arch3 <.Machine$double.eps,eps,arch3))%>%
    mutate(arch4=ifelse(arch4 <.Machine$double.eps,eps,arch4))
  return(MarkersCellsTMENS.ct%>%dplyr::select(c(arch1,arch2,arch3,arch4,value)))
  
} 

init <- c(0.1,0.1,0.1,0.1,1)
print("check2")
#### FOR MARKER CD20 IN B CELLS
x_cm<- get_tmens_coord_marker(MarkersCellsTMENs,"B","CD20",eps)
#x_cm2 <- get_tmens_coord_marker("CD4-T","PD1")
print("Starting Gradient descent 1")
start_time <- Sys.time()
estim_params <- optim(par = init, fn = func,data = sample_n(x_cm,1000),method = "L-BFGS-B",lower = c(-100,-100,-100,-100,0),upper = c(rep(100,NB_TMENs),Inf),
                      control = list(trace = 3,
                                     maxit = 1000,
                                     ndeps = rep(1e-4,NB_TMENs+1)))
write.csv(estim_params$par,"./fittedParamsGammaKm_CD20_B_reduced.csv")
end_time <- Sys.time()
print('output 1 written')
print(paste0("TIME 1: ",end_time-start_time))

#### FOR MARKER PD1 IN CD8-T CELLS
x_cm2<- get_tmens_coord_marker(MarkersCellsTMENs,"CD8-T","PD1",eps)
#x_cm2 <- get_tmens_coord_marker("CD4-T","PD1")
print("Starting Gradient descent 2")
start_time <- Sys.time()
estim_params2 <- optim(par = init, fn = func,data = sample_n(x_cm2,1000),method = "L-BFGS-B",lower = c(-100,-100,-100,-100,0),upper = c(rep(100,NB_TMENs),Inf),
                      control = list(trace = 3,
                                     maxit = 1000,
                                     ndeps = rep(1e-4,NB_TMENs+1)))
write.csv(estim_params2$par,"./fittedParamsGammaKm_PD1_CD8T_reduced.csv")
end_time <- Sys.time()
print('output 2  written')
print(paste0("TIME 2: ",end_time - start_time))

