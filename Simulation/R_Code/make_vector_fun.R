temp<-c()

make_vector<-function(x1,x2,vect1,vect2,n_sites){
  
   output<-rep(x1,x2)
    
  for (i in 2:n_sites){
    temp<-rep(vect1[i],vect2[i])
    output<-append(output,temp)
    temp<-c()
  }
  return(output)
}
