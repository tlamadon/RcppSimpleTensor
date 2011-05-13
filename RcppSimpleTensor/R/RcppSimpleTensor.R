# TODO: add support for recursive formulations
# for simulating time series for example

# I create a multidimensional array
library(Rcpp)
library(inline)


# this function takes a term and 
# grabs informations about the tensor expression
# it returns the list of Tensors, the list of indexes
# and the list of scalars
RcppSimpleTensorGetArgs <- function(a,r) {

  # if we find a scalar, we append it
  if (length(a)==1) {

    # check if we have a scalar, ow it's just a number
    if (! is.numeric(a))
      r$S = unique(c(r$S,paste(a)));
    
    r$E = paste(a)
    return(r)
  } 

  # if we find a tensor we do:
  # 1) we check if we know about all the indices and append the new ones
  # 2) we rewrite the bracket expression using C pointers
  # 3) we return that string
  ind = ""
  if (length(a)>=3 & a[[1]] == '[') {
    
    tensorName = paste(a[[2]]) 

    # running through the indices
    for (i in c((length(a)-2):1)){

      index  = paste(a[[i+2]])
      indexn = paste(index,"_n",sep="")
      indexi = paste(index,"_i",sep="")

      # appending the indice if it's new
      r$I = unique(c(r$I,index)) 
        
      # storing that this indices is this dimension for this tensor
      r$D = rbind(r$D , data.frame(I=index, M=tensorName, dim = i))

      # using the size and index, create the Cpp pointer with correct value
      if (i == (length(a)-2)) 
        ind = indexi
      else
        ind = paste( indexi , " + " , indexn , "*( ", ind ,")",sep="")
    }

    # finally we can append the tensor name in the front
    STR = paste(tensorName, "[" ,ind, "]",sep="")

    r$E = STR
    r$M=unique(c(r$M,tensorName))
    return(r)
  }

  # if we encounter a binary operator
  if (length(a)==3) {

    op = paste(a[[1]])

    # process the first argument
    r = RcppSimpleTensorGetArgs( a[[2]],r)
    STR1 = r$E
    # process the second argument
    r = RcppSimpleTensorGetArgs( a[[3]],r)
    STR2 = r$E

    # compose the expression
    if (op %in% c("+","-","/","*",">","<","<=",">="))
      r$E  = paste( "(", STR1 ,")",op,"(",STR2,")",sep="");

    if (op %in% c("max","min"))
      r$E  = paste( "f" , op, "(", STR1 ,",",STR2,")",sep="");

    if (op ==  '^')
      r$E  = paste( "pow(", STR1 ,",",STR2,")",sep="");

    return(r)
  }

  if (length(a)==2) {
    r = RcppSimpleTensorGetArgs( a[[2]],r)
   
    if (a[[1]] == 'I') {
      r$E = paste("((", r$E , ")?1:0)",sep="")
    } 
    
    return(r)
  } 
}

# creates the c fucntion from expression
RcppSimpleTensor <- function(expr,verbose=FALSE) {

  # for testing
  # a = terms(R[i,j] ~ beta * B[i,k] * A[j,k])
  # a = terms(R[i,j] ~ I( beta > S[i,j]))
  # parse the expression usign terms
  a = terms(expr)

  if (a[[1]] != '~') {
    cat("the formula needs to be of the form R[...] ~ ...")
    return()
  }

  # extract right hand side
  r= list(D= data.frame(), S=c(), M=c() , I=c())
  RHS = RcppSimpleTensorGetArgs(a[[3]],r)
  LHS = RcppSimpleTensorGetArgs(a[[2]],r)

  # find indices to sum over
  indiceSum = setdiff(RHS$I,LHS$I)
  indiceOut = LHS$I

  # creating signature and
  # source head
  sig = c()
  sign = c()
  src=""
  src2=""
  srcloop=""

  # Adding the matrices/tensors
  for (i in RHS$M) {
    varn = i 
    sig  = c(sig , "numeric")
    sign = c(sign, paste(varn,"_",sep=""))
    src = paste(src , 'Rcpp::NumericVector ' , varn, '(' ,varn, '_);\r\n', sep="")
  }

  # Adding the scalars
  for (i in RHS$S) {
    varn = i
    sig  = c(sig , "double")
    sign = c(sign, paste(varn,"_",sep=""))
    src = paste(src , 'double ' , varn , ' = Rcpp::as<double> (', varn, '_);\r\n',sep="")
  }
 
  reduceSig = sig
  names(reduceSig) = sign;

  # Adding indices and sizes
  for (i in c(indiceOut,indiceSum)) {
    varn = i
    sig  = c(sig , "int")
    sign = c(sign, paste(varn,"_",sep=""))
    src = paste(src , 'int ' , varn , '_n = Rcpp::as<int> (', varn, '_);\r\n',sep="")
    src2 = paste(src2 , 'int ' , i , '_i;\r\n',sep="")
    srcloop = paste(srcloop,'for (',i,'_i=0; ',i,'_i<',i,'_n; ',i,'_i++)\r\n',sep="")
  }

  names(sig) = sign;

  # creating R sig
  Rsize=c()
  for (i in indiceOut) {
    Rsize = c(Rsize,paste(i,"_n",sep="")) 
  }
  Rsize = paste(Rsize,collapse="*")

  # pasting all together!
  CODE = src;
  CODE = paste(CODE,src2)
  CODE = paste(CODE,'Rcpp::NumericVector R(', Rsize, ');' ,"\r\n", sep="")
  CODE = paste(CODE,srcloop)
  CODE = paste(CODE,"{\r\n")
  CODE = paste(CODE, LHS$E , " = " , LHS$E,"+", sep="")
  CODE = paste(CODE, RHS$E, ";\r\n" , sep="")
  CODE = paste(CODE,"}\r\n")
  CODE = paste(CODE,"return R;\r\n")
  cat(CODE)

  tmpfun <- cxxfunction(sig, CODE, plugin="Rcpp",includes="#include <math.h>",verbose=verbose)

  # finally we wrap it into another function that will get the sizes 
  # automatically from the matrices

  # I construct again a function instead of using do.call which
  # seems very slow, I want to frontload as much as possible

  WRAPFUNC = paste("tmpFunWrap <- function(", paste(names(reduceSig),collapse=",") ,") {",";", sep="");
  # Adding automatic computation of the indices dimensions
  dd = RHS$D
  for (i in c(indiceOut,indiceSum)) {
    # find where to get the size from 
    ddt = dd[dd$I==i,]
    WRAPFUNC = paste(WRAPFUNC, i , "_ = dim(", ddt$M[1] ,"_)[", ddt$dim[1] ,"]" ,";",sep ="")  
  }
  # Calling the C function
  WRAPFUNC = paste(WRAPFUNC, "R = tmpfun(", paste(names(sig),collapse=",") ,");",sep ="")  
  # Reshape the answer
  dd = LHS$D
  WRAPFUNC = paste(WRAPFUNC, "dim(R) <- c(", paste(paste(dd$I[order(dd$dim)],"_",sep=""),collapse=","),") ;",sep="")
  WRAPFUNC = paste(WRAPFUNC, "return(R)}" ,sep="")

  if (verbose) {
    cat(WRAPFUNC)
  }

  eval(parse(text=WRAPFUNC))

  return(tmpFunWrap)
}

# myfun <- tensorFunc(R[i,j] ~ beta * B[i,j])
# for (i in c(1:300)) {
# R = myfun(A1,0.6)
# }


#plug<- function(){
#  settings<- getPlugin( "Rcpp" )
#  settings$env[["PKG_CPPFLAGS"]]<- " -O3 "

#}
#registerPlugin( "omp", plug )





