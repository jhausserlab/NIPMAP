rm(list=ls())
library(tidyverse)

L = 800 #um
nNiches = 4

colScheme =
  matrix(paste0("0x", c("46", "cb", "ec",
                        "00", "00", "00",
                        "ff", "00", "00",
                        "ff", "00", "df")), nrow=4, ncol=3, byrow=T) %>% 
  strtoi %>% matrix(nrow=4) / 255

########################################

# diffusion2D <- function(t, Y, par) {
#   y    <- matrix(nrow = nx, ncol = ny, data = Y)  # vector to 2-D matrix
#   dY   <- beta*(1-y)*y       # consumption
#   
#   ## diffusion in X-direction; boundaries=imposed concentration
#   Flux <- -D * rbind(y[1,]-y[nx,], (y[2:nx,] - y[1:(nx-1),]), y[1,]-y[nx,])/dx
#   dY   <- dY - (Flux[2:(nx+1),] - Flux[1:nx,])/dx
#   
#   ## diffusion in Y-direction
#   Flux <- -D * cbind(y[,1]-y[,ny], (y[,2:ny]-y[,1:(ny-1)]), y[,1]-y[,ny])/dy
#   dY   <- dY - (Flux[,2:(ny+1)] - Flux[,1:ny])/dy
#   
#   return(list(as.vector(dY)))
# }
# 
# ## parameters
# dy <- dx <- 1  # grid size
# D <- 1  # diffusion coeff, X- and Y-direction
# beta  <- 0.05     # consumption rate
# 
# nx <- 50
# ny <- nx
# y  <- matrix(runif(nx*ny), nrow = nx, ncol = ny)
# 
# ## model most efficiently solved with lsodes - need to specify lrw
# 
# library(deSolve)
# print(system.time(
#   ST3 <- ode.2D(y, times = 1:100, func = diffusion2D, parms = NULL,
#                 dimens = c(nx, ny), verbose = TRUE, names = "Y",
#                 lrw = 400000, atol = 1e-10, rtol = 1e-10, cyclicBnd = 2)
# ))
# 
# # summary of 2-D variable
# summary(ST3)
# 
# # plot output at t = 10
# t10 <-  matrix (nrow = nx, ncol = ny, 
#                 data = subset(ST3, select = "Y", subset = (time == 10)))
# 
# persp(t10, theta = 30, border = NA, phi = 70, 
#       col = "lightblue", shade = 0.5, box = FALSE)
# 
# zlim <- range(ST3[, -1])
# for (i in 2:nrow(ST3)) {
#   y <- matrix(nr = nx, nc = ny, data = ST3[i, -1])
#   filled.contour(y, zlim = zlim, main = i)
# }

#############################################################

diffusion2D <- function(t, Y, par) {
  y    <- array(Y, dim=c(nx, ny, nNiches))  # vector to 2-D matrix
  
  totY = apply(y, c(1,2), sum)
  totYrep = array(NA, c(nx, ny, nNiches))
  for (i in 1:nNiches) {
    totYrep[,,i] = totY
  }
#  dY   <- beta*(1-totYrep)*y^2       # reaction
#  dY   <- beta*(1-totYrep)*(y^2/(y^2+(1/nNiches)^2))       # reaction
  #dY   <- beta*(1-totYrep)*y       # reaction
  
  K = 1/nNiches
  dY   <- beta*y*(y^4/(y^5+K^5) - 1/(2*K)) #reaction
  
  ## diffusion in X-direction; boundaries=imposed concentration
  FluxX = array(NA, dim=c(nx+1, ny, nNiches))
  # FluxX[1,] <- -D * rbind(y[1,,]-y[nx,,], (y[2:nx,,] - y[1:(nx-1),,]), y[1,,]-y[nx,,])/dx
  FluxX[1,,] = y[1,,]-y[nx,,]
  FluxX[2:nx,,] = y[2:nx,,] - y[1:(nx-1),,]
  FluxX[nx+1,,] = y[1,,]-y[nx,,]
  FluxX = -D * FluxX / dx
  dY   <- dY - (FluxX[2:(nx+1),,] - FluxX[1:nx,,]) / dx
  
  ## diffusion in Y-direction
  # Flux <- -D * cbind(y[,1,]-y[,ny,], (y[,2:ny,]-y[,1:(ny-1),]), y[,1,]-y[,ny,])/dy
  FluxY = array(NA, dim=c(nx, ny+1, nNiches))
  FluxY[,1,] = y[,1,]-y[,ny,]
  FluxY[,2:ny,] = y[,2:ny,] - y[,1:(ny-1),]
  FluxY[,ny+1,] = y[,1,]-y[,ny,]
  FluxY = -D * FluxY / dy
  
  dY   <- dY - (FluxY[,2:(ny+1),] - FluxY[,1:ny,])/dy
  
  return(list(as.vector(dY)))
}

## parameters
# nx <- 50
# ny <- nx

#dx <- L / nx   # grid size
dx = 16 #um
dy = dx

nx = L / dx
ny = nx

D <- .01 * L  # diffusion coeff, X- and Y-direction
D <- .045 * L  # diffusion coeff, X- and Y-direction
D <- .05 * L  # diffusion coeff, X- and Y-direction
D = 40
# D <- 0 * L  # diffusion coeff, X- and Y-direction
beta  <- 1     # max proliferation rate
fracLastNiche = .25 #abundance of last niche

y <- array(runif(nx*ny*nNiches), dim=c(nx,ny,nNiches))
y[,,nNiches] = y[,,nNiches] * (nNiches * fracLastNiche)^(1/(nNiches-1))
## model most efficiently solved with lsodes - need to specify lrw
mean(y[,,4] > apply(y[,,1:3], 1:2, max))

alphaInit = apply(y, c(1:2), function(x) { x == max(x) }) %>% aperm(c(2,3,1))
dim(alphaInit)
mean(alphaInit[,,4])
alphaInit = as.numeric(alphaInit)
#alphaInit = alphaInit * runif(length(alphaInit))

library(deSolve)
print(system.time(
  ST3 <- ode.2D(alphaInit, times = 1:1e2, func = diffusion2D, parms = NULL,
                nspec=nNiches, dimens = c(nx, ny), verbose = TRUE, names = paste0("niche", 1:nNiches),
                lrw = 500e6, atol = 1e-10, rtol = 1e-10, cyclicBnd = c(1,2))
))

# summary of 2-D variable
# summary(ST3)
par(mfrow=c(2,2))
plot(sapply(2:nrow(ST3), function(i) { mean(sqrt((ST3[i,] - ST3[i-1,])^2)) }), ylab="Cauchy diff")
plot(apply(ST3, 1, function(x) {array(x, dim=c(nx, ny, nNiches))[,,nNiches] %>% sd }), ylab="sd")
plot(apply(ST3, 1, function(x) {array(x, dim=c(nx, ny, nNiches))[,,nNiches] %>% mean }), ylab="mean", ylim=c(0,.6))
plot(apply(ST3, 1, function(x) {array(x, dim=c(nx, ny, nNiches))[,,nNiches] %>% max }), ylab="max", ylim=c(0,.6))

# # plot output at ss
# t10 <-  matrix (nrow = nx, ncol = ny, 
#                 data = subset(ST3, select = "niche4", subset = (time == nrow(ST3))))
# persp(t10, theta = 30, border = NA, phi = 70, 
#       col = "lightblue", shade = 0.5, box = FALSE)
# summary(t10 %>% as.numeric)
# 
# zlim <- c(0,1) #range(ST3[, -1])
# par(mfrow=c(3,3))
#  for (i in round(seq(1, nrow(ST3), len=3*3))) {
#    t10 <-  matrix (nrow = nx, ncol = ny, 
#                    data = subset(ST3, select = "niche4", subset = (time == i)))
#    persp(t10, theta = 30, border = NA, phi = 70, 
#          col = "lightblue", shade = 0.5, box = FALSE)
#  }

endState = array(ST3[nrow(ST3),], dim=c(nx, ny, nNiches))
apply(endState, 3, mean)
apply(endState, 3, max)
image(endState[,,4])

stateToRGB = function(endState, colScheme) {
  ix = 1
  iy = 1
  imgL = map(1:nx, function(ix) {
    map(1:ny, function(iy) {
      w = matrix(endState[ix,iy,] / sum(endState[ix,iy,]), nrow=1)
      (w %*% colScheme) %>% as.numeric()
    })
  })
  
  img = array(as.numeric(unlist(imgL)), dim=c(3, ny, nx))
  img = aperm(img, c(3,2,1))
}

array(ST3[1,], dim=c(nx, ny, nNiches))

for (i in 1:9) {
  endState = array(ST3[i,], dim=c(nx, ny, nNiches))
  print(endState[1:5,1:5,nNiches])
}

par(mfrow=c(3,3))
for (i in round(seq(1, nrow(ST3), len=3*3))) {
#for (i in 1:9) {
  endState = array(ST3[i,], dim=c(nx, ny, nNiches))
  endState[1:3,1:3,]
  rgbImg = stateToRGB(endState, colScheme) 
  rgbImg[rgbImg>1] = 1
  plot.new()
  plot.window(xlim=c(0,1),ylim=c(0,1))
  rasterImage(rgbImg, 0,0,1,1, interpolate=T) 
}

mean(endState[,,4] > apply(endState[,,1:3], 1:2, max))

# for (i in 2:nrow(ST3)) {
#   y <- matrix(nr = nx, nc = ny, data = ST3[i, -1])
#   filled.contour(y, zlim = zlim, main = i)
# }

### Calibrate fracLastNiche to generate alpha_i maps with a range of TLS abundance

library(furrr)
plan(multisession, workers = 7)

# This does run in parallel!
# future_map(c("hello", "world"), ~.x)

fracLastNiches = exp(seq(log(.04), log(.25), len=10))

setwd("/home/jean/projects/lab/manuscripts/NIPMAP")
resTib = future_map(fracLastNiches, function(fracLastNiche) {
  # fracLastNiche = .25 #abundance of last niche
  
  y <- array(runif(nx*ny*nNiches), dim=c(nx,ny,nNiches))
  y[,,nNiches] = y[,,nNiches] * (nNiches * fracLastNiche)^(1/(nNiches-1))
  ## model most efficiently solved with lsodes - need to specify lrw
  
  alphaInit = apply(y, c(1:2), function(x) { x == max(x) }) %>% aperm(c(2,3,1))
  alphaInit = as.numeric(alphaInit)
  
  ST3 <- ode.2D(alphaInit, times = 1:1e2, func = diffusion2D, parms = NULL,
                nspec=nNiches, dimens = c(nx, ny), verbose = TRUE, names = paste0("niche", 1:nNiches),
                lrw = 1e6, atol = 1e-10, rtol = 1e-10, cyclicBnd = c(1,2))
  endState = array(ST3[nrow(ST3),], dim=c(nx, ny, nNiches))
  nicheArea = mean(endState[,,nNiches] > apply(endState[,,1:(nNiches-1)], 1:2, max))
  write_rds(endState, file=sprintf("alphaMap_area%.4f_len%.1f.rds", nicheArea, L))
  return(tibble(fracNiche=fracLastNiche, nicheArea=nicheArea))
}) %>% bind_rows()
#write_rds(resTib, file="resTib.rds")
ggplot(resTib) + geom_line(aes(x=fracNiche, y=nicheArea)) + scale_x_log10() + scale_y_log10()

## Now sample cells 100 times according to each alpha, so we have data from 1..100 MIBI samples at each niche abundance
nSites = round(L^2 / (pi*25^2))

sitesT = tibble(x=runif(nSites)*L, y=runif(nSites)*L)
# select sites
# find alpha at the site
# from alpha and B (ask Anissa), determine nCells (do some Gaussian kernel maths) at the site and their type 
# -> cell type counts per site

# repeat for 1..100 samples and 10 different niche prevalances
exp(seq(log(1), log(1000), len=10))

 -> 10 niche prevalances x 10 dataset sizes x 10 repeats (re-sample cell count matrix per site)
 - for each one, ask Anissa if she finds a niche within eps of the TLS niche

# send to Anissa 
