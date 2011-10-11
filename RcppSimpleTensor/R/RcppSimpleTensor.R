# TODO: add support for recursive formulations
# for simulating time series for example

# TODO: use of tensors inside formula (i,j)
# TODO: return a scalar
# TODO: use the shape of a tensor supplied as extra arguments for indexes
# TODO: create inline tensor function
# notation could be : 
#   TI(expression ,  k+l|i+j       ) or
#   TI(expression ,  shape=G[i,j,k])


# I create a multidimensional array
require(Rcpp)
require(inline)
require(digest)

# check out ?methods to overload $


# this function takes a term and 
# grabs informations about the tensor expression
# it returns the list of Tensors, the list of indexes
# and the list of scalars
RcppSimpleTensorGetArgs <- function(a,r) {

  # if we find a scalar, we append it
  # unless it's in i,j,k,l, in which case we replace it with i_i ot j_i
  if (length(a)==1) {

    # check if we have a scalar, ow it's just a number
    vv = paste(a)
    if ( vv %in% c('i','j','k','l')) {
      vv= paste(vv,'_i',sep='')
    } else if (!is.numeric(a)) {
      r$S = unique(c(r$S,paste(a)))
    }

    r$E = vv 
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
      at = a[[i+2]]

      # grab a simple lag
      if (length(at)>1) {
        indexlag = as.numeric(at[[3]])
        index  = paste(at[[2]])      
      } else {
	indexlag=0
        index  = paste(at)      
      }

      indexn = paste(index,"_n",sep="")
      indexi = paste(index,"_i",sep="")

      # appending the indice if it's new
      r$I = unique(c(r$I,index)) 
        
      # storing that this indices is this dimension for this tensor
      r$D = rbind(r$D , data.frame(I=index, M=tensorName, dim = i, lag=indexlag))

      # create lag string
      if (indexlag>0) {
        LAG_STR = paste('-',indexlag)
      } else {
	LAG_STR = ''
      }

      # using the size and index, create the Cpp pointer with correct value
      if (i == (length(a)-2)) 
        ind = paste(indexi,LAG_STR)
      else
        ind = paste( indexi ,LAG_STR, " + " , indexn , "*( ", ind ,")",sep="")
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
    if (op %in% c("+","-","/","*",">","<","<=",">=",'||','&&'))
      r$E  = paste(STR1,op,STR2,sep=" ");
      #r$E  = paste( "(", STR1 ,")",op,"(",STR2,")",sep="");

    if (op %in% c("max","min"))
      r$E  = paste( "f" , op, "(", STR1 ,",",STR2,")",sep="");

    if (op %in%  c('^','pow'))
      r$E  = paste( "pow(", STR1 ,",",STR2,")",sep="");

    return(r)
  }

  # we found a unary operator
  if (length(a)==2) {
    r = RcppSimpleTensorGetArgs( a[[2]],r)
   
    if (a[[1]] == 'I') {
      r$E = paste("((", r$E , ")?1.0:0.0)",sep="")
    }

    if (paste(a[[1]]) %in% c('exp')) {
      r$E = paste(a[[1]],"(", r$E , ")",sep="")
    }

    if (paste(a[[1]]) %in% c('sign')) {
      r$E = paste("copysignf(1.0,", r$E , ")",sep="")
    }

    if (paste(a[[1]]) %in% c('abs')) {
      r$E = paste('f',a[[1]],"(", r$E , ")",sep="")
    }

    if (a[[1]] == '(') {
      r$E = paste("(", r$E , ")",sep="")
    }
    
    return(r)
  } 
}

# creates the c fucntion from expression
RcppSimpleTensor <- function(expr,cache=TRUE,verbose=FALSE,struct=FALSE) {

  # look if we already have this compiled localy
  fhash = digest(paste(expr,collapse=""))
  filename =paste('rcpptensor',fhash,sep="")

  # for testing
  # a = terms(R[i,j] ~ beta * B[i,k] * A[j,k])
  # a = terms(R[i,j] ~ I( beta > S[i,j]))
  # parse the expression usign terms
  if (is.character(expr)) {
    a = terms(as.formula(expr))
  } else {
    a = terms(expr)
  }

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
    if (varn == 'R') { next} #skip R
    sig  = c(sig , "numeric")
    sign = c(sign, paste(varn,"_",sep=""))
    src = paste(src , 'Rcpp::NumericVector ' , varn, '(' ,varn, '_);\n', sep="")
  }

  # Adding the scalars
  for (i in RHS$S) {
    varn = i
    sig  = c(sig , "double")
    sign = c(sign, paste(varn,"_",sep=""))
    src = paste(src , 'double ' , varn , ' = Rcpp::as<double> (', varn, '_);\n',sep="")
  }
 
  reduceSig = sig
  names(reduceSig) = sign;
  varframe= RHS$D

  if (verbose) {
    print(varframe)
  }

  # Adding indices and sizes
  for (i in c(indiceOut,indiceSum)) {
    varn = i
    # get highest lag
    lag = max(varframe$lag[varframe$I==varn])

    sig  = c(sig , "int")
    sign = c(sign, paste(varn,"_",sep=""))
    src = paste(src , 'int ' , varn , '_n = Rcpp::as<int> (', varn, '_);\n',sep="")
    src2 = paste(src2 , 'int ' , i , '_i;\n',sep="")
    srcloop = paste(srcloop,'for (',i,'_i=',lag,'; ',i,'_i<',i,'_n; ',i,'_i++)\n',sep="")
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
  CODE = paste(CODE,'Rcpp::NumericVector R(', Rsize, ');' ,"\n", sep="")
  CODE = paste(CODE,srcloop)
  CODE = paste(CODE,"{\n")
  CODE = paste(CODE, LHS$E , " = " , LHS$E,"+ \\\n", sep="")
  
  # we need to break the line, let's do it every 60 charaters
  # but we need to break on a blank space!
  LINE_LENGTH= 60
  LEFT_STR = strsplit(RHS$E," ")
  lcount = 0
  for (s in LEFT_STR[[1]]) {

    if (lcount>LINE_LENGTH) {
      CODE = paste(CODE, "\\\n" , sep="")
      lcount =0
    }

    CODE = paste(CODE, s , sep="")
    lcount = lcount + nchar(s)
  }
  CODE = paste(CODE, "; \n" , sep="")
  
  #CODE = paste(CODE, RHS$E, ";\n" , sep="")

  CODE = paste(CODE,"}\n")
  CODE = paste(CODE,"return R;\n")
  if (verbose) {
    #cat(CODE)
  }

  tmpfun <- mycxxfunction(sig, CODE,plugin="Rcpp",file=filename,includes="#include <math.h>",verbose=verbose,cache=cache)

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
  
  if (struct) {
    res = list()
    res$inFunc = tmpfun
    res$sig    = sig
    res$wrapFunc = tmpFunWrap
    res$RHS = RHS
    return(res)
  }else {
    return(tmpFunWrap)
  }

}

#}
#registerPlugin( "omp", plug )
mycxxfunction <- function (sig = character(), body = character(), plugin = "default", 
    includes = "", settings = getPlugin(plugin), file= basename(tempfile()) , ..., verbose = FALSE , cache = FALSE) 
{
    #f <- basename(tempfile())
    f = file
    if (!is.list(sig)) {
        sig <- list(sig)
        names(sig) <- f
        if (!length(body)) 
            body <- ""
        names(body) <- f
    }
    if (length(sig) != length(body)) 
        stop("mismatch between the number of functions declared in 'sig' and the number of function bodies provided in 'body'")
    signature <- lapply(sig, function(x) {
        if (!length(x)) {
            ""
        }
        else {
            paste(sprintf("SEXP %s", names(x)), collapse = ", ")
        }
    })
    decl <- lapply(1:length(sig), function(index) {
        sprintf("SEXP %s( %s) ;", names(signature)[index], signature[[index]])
    })
    def <- lapply(1:length(sig), function(index) {
        sprintf("\nSEXP %s( %s ){\n%s\n}\n", names(signature)[index], 
            signature[[index]], if (is.null(settings$body)) 
                body[[index]]
            else settings$body(body[[index]]))
    })
    settings_includes <- if (is.null(settings$includes)) 
        ""
    else paste(settings$includes, collapse = "\n")
    code <- sprintf("\n// includes from the plugin\n%s\n\n// user includes\n%s\n\n// declarations\nextern \"C\" {\n%s\n}\n\n// definition\n%s\n\n", 
        settings_includes, paste(includes, collapse = "\n"), 
        paste(decl, collapse = "\n"), paste(def, collapse = "\n"))
    if (!is.null(env <- settings$env)) {
        do.call(Sys.setenv, env)
        if (isTRUE(verbose)) {
            cat(" >> setting environment variables: \n")
            writeLines(sprintf("%s = %s", names(env), env))
        }
    }
    LinkingTo <- settings$LinkingTo
    if (!is.null(LinkingTo)) {
        paths <- .find.package(LinkingTo, quiet = TRUE)
        if (length(paths)) {
            flag <- paste(inline:::paste0("-I\"", paths, "/include\""), 
                collapse = " ")
            Sys.setenv(CLINK_CPPFLAGS = flag)
            if (isTRUE(verbose)) {
                cat(sprintf("\n >> LinkingTo : %s\n", paste(LinkingTo, 
                  collapse = ", ")))
                cat("CLINK_CPPFLAGS = ", flag, "\n\n")
            }
        }
    }
    if (isTRUE(verbose)) {
        writeLines(" >> Program source :\n")
        writeLines(inline:::addLineNumbers(code))
    }
    language <- "C++"

    # if tensor folder not here, we create it
    if (file.exists('.tensor')==FALSE){
        dir.create('.tensor') 
    }

    libLFile <- mycompileCode(f, code, language = language, verbose = verbose,dir  = paste(getwd(),'/.tensor/',sep=""),cache = cache)
    #cleanup <- function(env) {
    #    if (f %in% names(getLoadedDLLs())) 
    #        dyn.unload(libLFile)
    #    unlink(libLFile)
    #}
    #reg.finalizer(environment(), cleanup, onexit = TRUE)
    res <- vector("list", length(sig))
    names(res) <- names(sig)
    res <- new("CFuncList", res)

    #load only if not loaded yet
    dlllist = getLoadedDLLs()
    if (f %in%  names(dlllist)) {
      DLL <- dlllist[f]
      DLL = DLL[[1]]
    } else {
      DLL <- dyn.load(libLFile)
    }

    for (i in seq_along(sig)) {
        res[[i]] <- new("CFunc", code = code)
        fn <- function(arg) {
            NULL
        }
        args <- formals(fn)[rep(1, length(sig[[i]]))]
        names(args) <- names(sig[[i]])
        formals(fn) <- args
        body <- quote(.Call(EXTERNALNAME, ARG))[c(1:2, rep(3, 
            length(sig[[i]])))]
        for (j in seq(along = sig[[i]])) body[[j + 2]] <- as.name(names(sig[[i]])[j])
        body[[1L]] <- .Call
        body[[2L]] <- getNativeSymbolInfo(names(sig)[[i]], DLL)$address
        body(fn) <- body
        res[[i]]@.Data <- fn
    }
    rm(j)
    convention <- ".Call"
    if (identical(length(sig), 1L)) 
        res[[1L]]
    else res
}

mycompileCode <- function (f, code, language, verbose, dir = tmpdir(),cache=FALSE) 
{
    wd = getwd()
    on.exit(setwd(wd))
    if (.Platform$OS.type == "windows") {
        libCFile <- paste(dir, "/", f, ".EXT", sep = "")
        libLFile <- paste(dir, "/", f, ".dll", sep = "")
        libLFile2 <- paste(dir, "/", f, ".dll", sep = "")
    }
    else {
        libCFile <- paste(dir, "/", f, ".EXT", sep = "")
        libLFile <- paste(dir, "/", f, ".so", sep = "")
        libLFile2 <- paste(dir, "/", f, ".sl", sep = "")
    }

    if (cache & file.exists(libLFile))
      return(libLFile)

    extension <- switch(language, `C++` = ".cpp", C = ".c", Fortran = ".f", 
        F95 = ".f95", ObjectiveC = ".m", `ObjectiveC++` = ".mm")
    libCFile <- sub(".EXT$", extension, libCFile)
    write(code, libCFile)
    if (file.exists(libLFile)) 
        file.remove(libLFile)
    if (file.exists(libLFile2)) 
        file.remove(libLFile2)
    setwd(dirname(libCFile))
    errfile <- paste(basename(libCFile), ".err.txt", sep = "")
    cmd <- paste(R.home(component = "bin"), "/R CMD SHLIB ", 
        basename(libCFile), " 2> ", errfile, sep = "")
    if (verbose) 
        cat("Compilation argument:\n", cmd, "\n")
    compiled <- system(cmd, intern = !verbose)
    errmsg <- readLines(errfile)
    unlink(errfile)
    writeLines(errmsg)
    setwd(wd)
    if (!file.exists(libLFile) && file.exists(libLFile2)) 
        libLFile <- libLFile2
    if (!file.exists(libLFile)) {
        cat("\nERROR(s) during compilation: source code errors or compiler configuration errors!\n")
        cat("\nProgram source:\n")
        code <- strsplit(code, "\n")
        for (i in 1:length(code[[1]])) cat(format(i, width = 3), 
            ": ", code[[1]][i], "\n", sep = "")
        stop(paste("Compilation ERROR, function(s)/method(s) not created!", 
            paste(errmsg, collapse = "\n")))
    }
    return(libLFile)
}

#Sum <- function(expr,

TI <- function(arga,argb) {

  dims   = deparse(substitute(argb))
  dims   = gsub('\\+',',',dims)
  tsor   = deparse(substitute(arga))
  TENSOR = paste('R[',dims,'] ~ ',tsor,sep='',collapse='')
  
  rr = RcppSimpleTensor(TENSOR,struct=TRUE)

  # just need to return rr$wrapFunc applied to the arrays in the parent frame
  # we call funcWrap with the list of M and S from rr$RHS
  mycall = paste( 'rr$wrapFunc(',paste(
                               paste(c(rr$RHS$M,rr$RHS$S),sep='',collapse=","),
                               sep=','), ')')

  return(eval(parse(text=mycall)))
}

