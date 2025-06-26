#' @export
sparsityForThird = function(hessian, data, 
  Nbeta = nrow(data$XTp), Ngamma= nrow(data$ATp)) {

if(is.data.frame(hessian))
  hessian = Matrix::sparseMatrix(i=hessian$i, j=hessian$j, index1=FALSE)

if(nrow(hessian) <= Ngamma ) 
  warning("hessian must be from the full set of parameters and random effects")

Sgamma = seq(from=Nbeta+1, len=Ngamma)

hessianParametersGamma = as(hessian[Sgamma, -Sgamma], "TsparseMatrix")
# hessian of random effects, lower triangle only, column format
hessianRandom = as(
    as(hessian[Sgamma, Sgamma], "symmetricMatrix"),
  'CsparseMatrix')

# parameter-random cross hessian, triplet format
# to do: assum zero in tensor 3 if hessian zero for parameter-i combination
parametersGamma = cbind(
  gamma=hessianParametersGamma@i,
  gammaInParams = hessianParametersGamma@i + Nbeta, 
  param=hessianParametersGamma@j, 
  paramInParams=hessianParametersGamma@j + 
    (hessianParametersGamma@j>Nbeta)*(Ngamma),
  NinHessian = diff(hessianRandom@p)[1+hessianParametersGamma@i],
  PinHessian = hessianRandom@p[1+hessianParametersGamma@i]
  )
parametersGamma = cbind(parametersGamma, 
  PinThird = c(0, cumsum(parametersGamma[,'NinHessian'])[
    -nrow(parametersGamma)]))
storage.mode(parametersGamma) = 'integer'

# indices for the full third deriv sparse tensor
parametersGamma2 = data.table::as.data.table(parametersGamma)
parametersGammaBig <- parametersGamma2[
  , {
      rep_times   <- seq_len(NinHessian)        # 1…NinHessian
      newPin      <- PinHessian + rep_times - 1 # vector con incrementos
      .(PinHessian = newPin)                    # devuelve el bloque expandido
    },
  by = .(param, paramInParams, gamma, gammaInParams)
]

parametersGammaBig$i = 
  hessianRandom@i[1+parametersGammaBig$PinHessian]
parametersGammaBig$iInParams = parametersGammaBig$i + Nbeta


# quantities needs to subtarct off T_ijk from taylor coefficients
indexForDiag = as.matrix(data.table::melt(parametersGammaBig,
  id.vars = 'gammaInParams', value.name = 'index',
  measure.vars = c('paramInParams','iInParams'))[,c('gammaInParams','index')])
storage.mode(indexForDiag) = 'integer'
indexForDiag = indexForDiag[!duplicated(indexForDiag), ]
indexForDiag = Matrix::sparseMatrix(
  j=indexForDiag[,'gammaInParams'],
  i=indexForDiag[,'index'],
  index1=FALSE, dims = dim(hessian)
  )



  return(list(
    parametersGamma = as.data.frame(parametersGamma),
    sparsity = list(i = hessianRandom@i, p=hessianRandom@p),
    full = parametersGammaBig,
    indexForDiag = list(i=indexForDiag@i, p=indexForDiag@p)
    ))
}

