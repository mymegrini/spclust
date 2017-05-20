
random = function(w){
    n = length(w)
    # discrete cumulative distribution 
    p = cumsum(w)
    # simulate probability using Inverse Transform Sampling
    u = runif(1,0,p[n])
    # find the index of the inverse using dichotomic search
    a = 0
    b = n
    # loop invariant : a<b and p[a]<=u<p[b]
    while(b-a>1){
        c = a + (b-a)%/%2
        if (p[c]>u){
            b = c
        } else {
            a = c
        }
    }
    return(b)
}

dist2 = function(x, y){
    return (t(x-y)%*%(x-y))
}

iclusters = function(X, C){
    # sample size
    n = ncol(X)
    # cluster count
    k = ncol(C)
    # Cindex will contain the index of the nearest center
    Cindex = vector("double", n)
    # partitioning sample over cluster centers
    for (x in 1:n){
        d = dist2(X[,x],C[,1])
        Cindex[x] = 1
        for (c in 2:k){
            if (dist2(X[,x],C[,c]) < d){
                d = dist2(X[,x],C[,c])
                Cindex[x] = c
            }
        }
    }
    return(Cindex)
}

clusters = function(X, index){
    # cluster count
    k = ncol(C)
    # constructing the list of clusters
    clist = vector("list", k)
    for (c in 1:k){
        cindex = Cindex==c
        clist[[c]] = X[,cindex]
    }
    return(clist)
}

d2seed = function(X, k){
    # sample dimensions
    d = nrow(X)
    n = ncol(X)
    # allocate a dxk matrix for the centers
    C = matrix(NA_real_, d, k)
    # choose first center uniformely from the sample
    i = floor(runif(1,0,n))
    C[,1] = X[,i+1]
    # calculate initial weight distribution
    d2 = vector("double", n)
    for (i in 1:n){
        d2[i] = dist2(X[,i],C[,1])
    }
    # choose the next (k-1) remaining centers using the D² method
    for (c in 2:k){
        # select the next center using auxiliary function
        r = random(d2)
        C[,c] = X[,r]
        # update weight distribution
        for (i in 1:n){
            d2[i] = min(d2[i], dist2(X[,i],C[,c]))
        }
    }
    return(C)
}

kmpp = function(X, k, iter=0){
    # sample size
    n = ncol(X)
    # initialization of cluster centers
    C = d2seed(X,k)
    # optimization loop
    repeat{
        # Cindex will contain the index of the nearest scenter
        Cindex = vector("double",n)
        # partitioning sample over cluster centers
        for (x in 1:n){
            d = dist2(X[,x],C[,1])
            Cindex[x] = 1
            for (c in 2:k){
                if (dist2(X[,x],C[,c]) < d){
                    d = dist2(X[,x],C[,c])
                    Cindex[x] = c
                }
            }
        }
        # updating cluster centers
        convergence = TRUE
        for (c in 1:k){
            index = Cindex==c
            size = sum(index)
            cluster = X[,index]
            if (size>1){
                newcenter = rowSums(cluster)/size
            } else {
                newcenter = cluster
            }            
            if (all(newcenter!=C[,c])){
                convergence = FALSE
                C[,c] = newcenter
            }
        }
        iter = iter-1
        if(convergence || iter==0){
            break
        }
    }
    return(C)
}

egraph = function(d,e){
    n=nrow(d)
    for(i in 1:n){
        for(j in i:n){
            if(d[i,j]>e){
                d[i,j]=0
                d[j,i]=0
            }
        }
    }
    return(d)
}

knearest = function(k,d,mutual=FALSE){
    n=nrow(d)
    p=ncol(d)
    #Matrix n*p of zero
    M=matrix(data=numeric(n*p),ncol=p,nrow=n)
    for(i in 1:n){
        #Iterate the number of point that we have to connect
        kk=k
        while(kk!=0){
            #Cordinates of the nearest point of i
            #Start at (i,1)
            ii=i
            jj=1
            if(mutual){
                min2=d[1,i] 
                min=d[i,1] 
            }
            else{
                min=min(d[i,1],d[1,i])
            }
            #Found the cordinates of the nearest point of i
            for(j in 1:i){
                if(mutual){
                    if(d[i,j]<min && d[j,i]<min2){
                        min=d[i,j]
                        min2=d[j,i]
                        ii=i
                        jj=j
                    }
                }
                else{
                    if(d[i,j]<min || d[j,i]<min){
                        if(d[i,j]>d[j,i]){ d[i,j]=d[j,i] }
                        min=d[i,j]
                        ii=i
                        jj=j
                    }
                    
                }
            }
            kk=kk-1
            #M take the nearest point which is at (ii,jj)
            #M is symetric
            M[ii,jj]=M[jj,ii]=d[ii,jj]
        }
    }
    return(M)
}

degree = function(W){
    n=nrow(W)
    #Allocate a zero matrix length n*n
    D = matrix(data=numeric(n*n), nrow=n, ncol=n)
    for(i in 1:n){
        D[i,i]=sum(W[i,])
    }
    return(D)
}

laplacian = function(W,normalize=0){
    #length of the matrix W
    n=ncol(W)
    #degree matrix
    D=degree(W)
    L=D-W
    #Normalized Laplacian Lsym
    if(normalize==1){
        RD=solve(D^(1/2))
        return(RD*L*RD)
    }
    #Normalized Laplacian Lrw
    if(normalize==2){
        return(solve(D)*L)
    }
    #Unormalized Laplacian of W
    return(L)
}

eigenvectors = function(L,k){
    #Compute all the eigenvalues and eigenvectors of L
    X=eigen(L)
    #Take the first k eigenvectors
    t=length(X$values)
    return(X$vectors[,(t-k+1):t])
}

spec = function(S,k,normalize=0){
    # Weighted adjacency matrix of S
    W=S
    # Compute the normalized Laplacian Lrw
    L=laplacian(W,normalize)
    # Compute the first k-eigenvectors of Lrw
    U=eigenvectors(L,k)
    # Normalize the rows
    if(normalize==2){
        n = nrow(U)
        N=(rowSums(U^2)^(1/2))
        for(i in 1:n){
            if (N[i]>0){
                for(j in 1:k){
                    U[i,j]=U[i,j]/N[i]
                }
            }
        }
    }
    # Transpose the matrix U
    Y=t(U)
    # Cluster the points with k-means algorithm
    C=kmpp(Y,k)
    return(iclusters(Y,C))
}

pclusters = function(X, k, index, names, title, col = NULL, shape=1){
    # dataset dimension
    d = nrow(X)
    # computing boundaries for the data
    blist = vector("list",d)
    for(i in 1:d){
        blist[[i]] = c(min(X[i,]),max(X[i,]))
    }
    # constructing the list of clusters
    clist = vector("list", k)
    for (c in 1:k){
        bool = index==c
        clist[[c]] = X[,bool]
    }
    # plotting colors
    if (length(col)==0){
        col = c(2:(k+1))
    }
    # adjusting margins
    par(oma=c(2,2,5,2), mar=(c(0,0,0,0)+.5))
    # plotting in R²
    if (d==2){
        plot(clist[[1]][1,],clist[[1]][2,], col = col[1],
             pch=shape, xlab=names[1], ylab=names[2], 
             xlim = c(1,5), ylim = c(1,5))
        for (c in 2:k){
            points(clist[[c]][1,],clist[[c]][2,], col = col[c],
                   xlab=names[1], ylab=names[2], pch=shape)
        }
        title(title, outer=TRUE)
    }
    # plotting in higher dimensions
    else {
        par(mfrow=c(d,d))
        for (i in 1:d){
            for (j in 1:d){
                # plotting dimension name
                if (j==i){
                    plot(blist[[i]], blist[[i]], #ann = F,  
                         type = 'n', xaxt = 'n', yaxt = 'n',)
                    text(x = mean(blist[[i]]), y = mean(blist[[i]]),
                         paste(names[i]), cex = 1.25, col = "black")
                } else {
                    # plotting the clusters in [j,i] space
                    if(length(nrow(clist[[1]]))>0){
                        plot(clist[[1]][j,],clist[[1]][i,], ann=F,
                             col = col[1], xaxt='n', yaxt='n', 
                             pch=shape, xlim = blist[[j]], 
                             ylim = blist[[i]])
                        for (c in 2:k){
                            if(length(nrow(clist[[c]]))>0){
                                points(clist[[c]][j,],clist[[c]][i,], 
                                       ann=F, pch=shape, col = col[c],
                                       xaxt='n', yaxt='n')
                            }
                        }
                    } else {
                        plot(clist[[2]][j,],clist[[2]][i,], ann=F,
                         col = col[1], xaxt='n', yaxt='n', pch=shape, 
                         xlim = blist[[j]], ylim = blist[[i]])
                        for (c in 3:k){
                            if(length(nrow(clist[[c]]))>0){
                                points(clist[[c]][j,],clist[[c]][i,],
                                       ann=F, pch=shape, col = col[c],
                                       xaxt='n', yaxt='n')
                            }
                        }
                    }
                    
                    # drawing axes
                    if (i==d){
                        axis(1)
                    } else if (i==1){
                        axis(3)
                    }
                    if (j==d){
                        axis(4)
                    } else if (j==1){
                        axis(2)
                    }
                }
            }
        }
        title(title, outer=TRUE)
    }
}

# dataset : point cloud in R⁴
n = 100
X = c(runif(2*n,0,2),runif(n,2,4),runif(n,4,6))
dim(X) = c(4,n)
index = c(rep(1,50),rep(2,25),rep(3,25))

# clustering
kmcluster = iclusters(X,kmpp(X,3))
# plotting
pclusters(X,3,kmcluster,c("X1","X2","X3","X4"),
          "Kmeans Clustering : Simulated dataset")

# distance matrix
D = matrix(0,n,n)
for (i in 1:n){
    for (j in 1:n){
        D[i,j] = dist2(X[,i],X[,j])
    }
}
# ε-neighborhood graph
e = 3
S = egraph(D,e)
# clustering
spcluster = spec(S,3,0)
# plotting
pclusters(X,3,spcluster,c("X1","X2","X3","X4"),
          "Unnormalized Spectral Clustering : Simulated dataset")

d = read.csv(file="iris.data",head=F,sep=",")
d$V5<-factor(d$V5)
summary(d)

n = 150
X = t(matrix(c(d$V1,d$V2,d$V3,d$V4),nrow=n,ncol=4))
index = as.numeric(d$V5)
t = "Iris flower data (red=setosa,green=versicolor,blue=virginica)"
pclusters(X,3,index,c("Sepal length","Sepal width",
                      "Petal length","Petal width"),
          t ,shape=20)

# clustering
kmcluster = iclusters(X,kmpp(X,3))
# plotting
pclusters(X,3,kmcluster,c("Sepal length","Sepal width",
                          "Petal length","Petal width"),
          "Kmeans Clustering : Iris Flower data", shape=20)

# distance matrix
D = matrix(0,n,n)
for (i in 1:n){
    for (j in 1:n){
        D[i,j] = dist2(X[,i],X[,j])
    }
}
# ε-neighborhood graph
e = 1
S = egraph(D,e)
# clustering
spcluster = spec(S,3,0)
# plotting
pclusters(X,3,spcluster,c("Sepal length","Sepal width",
                          "Petal length","Petal width"), 
          "Unnormalized Spectral Clustering : Iris Flower data",
          shape=20)
