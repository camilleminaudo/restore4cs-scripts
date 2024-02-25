



  get_dxdy <- function(x,y){
    n <- length(y)
    dx <- (x[n]-x[1])/(n-1)
    dxdy <- NA*double(length = n)
    for (j in seq(2,n-1)){
      #dxdy[j] <- (y[j+1]-y[j-1])/(2*(dx))
      
      h1 = (x[j]-x[j-1]);
      h2 = (x[j+1]-x[j]);
      
      dxdy[j] = -(h2/(h1*(h1 +h2)))*y[j-1] - ((h1-h2)/(h1*h2))*y[j] + (h1/(h2*(h1 +h2)))*y[j+1]
      
    }
    # dxdy[1] <- (y[2]-y[1])/(x[2]-x[1])
    # dxdy[n] <- (y[n]-y[n-1])/(x[n]-x[n-1])
    
    
    
    h1 = (x[2]-x[1])
    h2 = (x[3]-x[2])
    
    dxdy[1] = -((2*h1 + h2)/(h1*(h1+h2)))*y[1] + ((h1 + h2)/(h1*h2))*y[2] - (h1/(h2*(h1 +h2)))*y[3]
    
    h1 = (x[n-1]-x[n-2])
    h2 = (x[n]-x[n-1])
    
    dxdy[n] = + (h2/(h1*(h1 +h2)))*y[n-2] - ((h1 + h2)/(h1*h2))*y[n-1] + ((h1 + 2*h2)/(h2*(h1+h2)))*y[n]
    
    
    return(dxdy)
  }

