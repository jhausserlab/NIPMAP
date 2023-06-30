
gramschmidt <- function(x) {
  x[x==0] <- 1e-100 # add a small number to avoid multiple with 0 -> otherwise we will have multiple NA in the result
  # Get the number of rows and columns of the matrix
  n <- ncol(x)
  m <- nrow(x)
  
  # Initialize the Q and R matrices
  q <- matrix(0, m, n)
  r <- matrix(0, n, n)
  
  for (j in 1:n) {
    v = x[,j] # Step 1 of the Gram-Schmidt process v1 = a1
    # Skip the first column
    if (j > 1) {
      for (i in 1:(j-1)) {
        r[i,j] <- t(q[,i]) %*% x[,j] # Find the inner product (noted to be q^T a earlier)
        # Subtract the projection from v which causes v to become perpendicular to all columns of Q
        v <- v - r[i,j] * q[,i] 
      }      
    }
    # Find the L2 norm of the jth diagonal of R
    r[j,j] <- sqrt(sum(v^2))
    # The orthogonalized result is found and stored in the ith column of Q.
    q[,j] <- v / r[j,j]
  }
  
  # Collect the Q and R matrices into a list and return
  qrcomp <- list('Q'=q, 'R'=r)
  return(qrcomp)
}


proj <- function(u, v) {
  return(as.vector((u%*% v)/(v%*%v))*v)
}

Faceprojection <- function(arche_list, niche_pccompo, niche_alfas,cols){
  
  combiniche <- combn(arche_list, 3) # generate all combinations of niches
  plist <- list()
  cols = cbind(color=cols, arche=c(arche_list)) %>% as.data.frame() #cbind(color=c("#FF00DF", "#0000DF", "#46CBEC", "#FF0000", "#000000"), arche=c(arche_list)) %>% as.data.frame()
  for (num in seq(1, ncol(combiniche))) {
    
    #print(paste("Round",num))
    niche_chosen <- combiniche[,num]
    col_chosen <- cols %>% filter(arche %in% c(niche_chosen))
    col_chosen <- as.vector(col_chosen$color)
    # print(niche_chosen)
    
    X <- niche_pccompo[niche_chosen,]%>%as.matrix()
    i = X[1,] %>% as.matrix()
    j = X[2,] %>% as.matrix()
    k = X[3,] %>% as.matrix()
    i_j = j-i
    i_k = k-i
    #X_new <- 
    i_jn = i_j / sqrt(sum(i_j^2))
    proj_ik_on_ij = sum(i_k*i_jn) * i_jn
    i_kOrtho = i_k - proj_ik_on_ij
    i_kOrthon = i_kOrtho / sqrt(sum(i_kOrtho^2))
    orthonbasis = cbind(i_jn, i_kOrthon)
    sum(orthonbasis[,1] * orthonbasis[,2])
    sum(orthonbasis[,1] * orthonbasis[,1])
    sum(orthonbasis[,2] * orthonbasis[,2])
    #gs <- gramschmidt(X_new %>% t()) # gs$Q is the orthonormal basis of the face
    #print("End gramschmidt computation")
    
    ijk_2D <- t(orthonbasis) %*% cbind(rep(0, ncol(X)), cbind(i_j, i_k))
    colnames(ijk_2D) = c("i", "j", "k")
    rownames(ijk_2D) = c("x", "y")
    
    # 
    # jk_2D <- proj(X_new, gs$Q) %>%
    #   rbind(i = rep(0, 2)) %>%
    #   `colnames<-`(c("x","y")) %>%
    #   `rownames<-`(c("ij_2D","ik_2D","i_2D")) %>%
    #   t() %>% as.data.frame() %>%
    #   select(i_2D,ij_2D,ik_2D) %>% t()%>%
    #   `rownames<-`(c(niche_chosen))
    
    #print("Begin compute alpha^hat")
    isCloseToFace = (apply(niche_alfas[,niche_chosen], 1, sum)>0.5)
    isCloseToFace[is.na(isCloseToFace)] = F
    mean(isCloseToFace)
    
    site_pos <- as.matrix(niche_alfas[isCloseToFace,niche_chosen]) %*% t(ijk_2D) %>% as.data.frame()
    # print("End compute the position of the projected site")
    
    # create scatter plot
    tijk_2D = as_tibble(ijk_2D %>% t)
    tijk_2D[,"arch"] = str_replace_all(niche_chosen,"arch","")
    
    p <- ggplot(site_pos %>%
                  rbind(t(ijk_2D)) %>% 
                  rownames_to_column("site_name"), aes(x=x, y=y)) +
      geom_point(alpha = 0.5, size = 0.01) + 
      # layer(geom = "point", 
      #       params = list(color="orange",size=5), 
      #       data=jk_2D%>%as.data.frame(), 
      #       stat = "identity", 
      #       position = "identity")+
      #annotate("text", x = c(jk_2D[,"x"]), y = c(jk_2D[,"y"]), label = c(niche_chosen) , color="orange", size=8, vjust = "inward", hjust = "inward") + 
      geom_point(data = tijk_2D%>%as.data.frame(), aes(x=x, y=y), size=5, color=col_chosen) +
      #geom_label(data = tijk_2D, aes(x=x, y=y, label=arch)) +
      ggtitle(paste(niche_chosen[1],niche_chosen[2],niche_chosen[3])) +
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank())
    #    xlab(niche_chosen[2]) + ylab(niche_chosen[3])
    
    plist[[num]] <- p
    
    # print("End round", num)
  }
  
  return(plist)
}
