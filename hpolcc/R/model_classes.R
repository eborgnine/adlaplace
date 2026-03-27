#' Model Classes for hpolcc
#'
#' This file defines S4 classes for different model types used in the hpolcc package.
#' These classes provide a structured way to handle model specifications and ensure
#' type safety when working with different model types.
#'
#' @name model_classes
NULL

#' @rdname model_classes
#' @export
setClass("iwp",
  representation = representation(
    var = "character",
    name = "character",
    f = "formula",
    p = "numeric",
    ref_value = "numeric",
    knots = "ANY",
    range = "ANY",
    init = "numeric",
    lower = "numeric",
    upper = "numeric",
    parscale = "numeric",
    random = "logical"
  ),
  prototype = list(
    var = character(0),
    name = character(0),
    f = formula(),
    p = numeric(0),
    ref_value = numeric(0),
    knots = NULL,
    range = NULL,
    init = numeric(0),
    lower = numeric(0),
    upper = numeric(0),
    parscale = numeric(0),
    random = TRUE
  )
)

#' @rdname model_classes
#' @export
setClass("hiwp",
         representation = representation(
           var = "character",
           name = "character",
           f = "formula",
           p = "numeric",
           ref_value = "numeric",
           knots = "ANY",
           range = "ANY",
           group_var = "character",
           groups = "integer",
           groups_string = "character",
           init = "numeric",
           lower = "numeric",
           upper = "numeric",
           parscale = "numeric",
           random = "logical"
         ),
         prototype = list(
           var = character(0),
           name = character(0),
           f = formula(),
           p = numeric(0),
           ref_value = numeric(0),
           knots = NULL,
           range = NULL,
           group_var = character(0),
           groups = integer(0),
           groups_string = character(0),
           init = numeric(0),
           lower = numeric(0),
           upper = numeric(0),
           parscale = numeric(0),
           random = TRUE
         )
)

#' @rdname model_classes
#' @export
setClass("linear",
  representation = representation(
    var = "character",
    name = "character",
    f = "formula",
    random = "logical"
  ),
  prototype = list(
    var = character(0),
    name = character(0),
    f = formula(),
    random = FALSE
  )
)

#' @rdname model_classes
#' @export
setClass("fpoly",
         representation = representation(
           var = "character",
           name = "character",
           f = "formula",
           p = "numeric",
           ref_value = "numeric",
           init = "numeric",
           lower = "numeric",
           upper = "numeric",
           parscale = "numeric",
           random = "logical"
         ),
         prototype = list(
           var = character(0),
           name = character(0),
           f = formula(),
           p = numeric(0),
           ref_value = numeric(0),
           init = numeric(0),
           lower = numeric(0),
           upper = numeric(0),
           parscale = numeric(0),
           random = FALSE
         )
)

#' @rdname model_classes
#' @export
setClass("rpoly",
         representation = representation(
           var = "character",
           name = "character",
           f = "formula",
           p = "numeric",
           ref_value = "numeric",
           sd = "numeric",
           random = "logical"
         ),
         prototype = list(
           var = character(0),
           name = character(0),
           f = formula(),
           p = numeric(0),
           ref_value = numeric(0),
           sd = numeric(0),
           random = TRUE
         )
)

#' @rdname model_classes
#' @export
setClass("hrpoly",
         representation = representation(
           var = "character",
           name = "character",
           f = "formula",
           p = "numeric",
           ref_value = "numeric",
           group_var = "character",
           groups = "integer",
           groups_string = "character",
           init = "numeric",
           lower = "numeric",
           upper = "numeric",
           parscale = "numeric",
           random = "logical"
         ),
         prototype = list(
           var = character(0),
           name = character(0),
           f = formula(),
           p = numeric(0),
           ref_value = numeric(0),
           group_var = character(0),
           groups = integer(0),
           groups_string = character(0),
           init = numeric(0),
           lower = numeric(0),
           upper = numeric(0),
           parscale = numeric(0),
           random = TRUE
         )
)

#' @rdname model_classes
#' @export
setClass("iid",
         representation = representation(
           var = "character",
           name = "character",
           f = "formula",
           init = "numeric",
           lower = "numeric",
           upper = "numeric",
           parscale = "numeric",
           random = "logical"
         ),
         prototype = list(
           var = character(0),
           name = character(0),
           f = formula(),
           init = numeric(0),
           lower = numeric(0),
           upper = numeric(0),
           parscale = numeric(0),
           random = TRUE
         )
)

#' Coerce hiwp to iwp
#'
#' @param from The hiwp object to coerce
#' @return An iwp object
#' @export
setAs("hiwp", "iwp",
       function(from) {
         new("iwp",
             var = from@var,
             name = from@name,
             f = from@f,
             p = from@p,
             ref_value = from@ref_value,
             knots = from@knots,
             range = from@range,
             init = from@init,
             lower = from@lower,
             upper = from@upper,
             parscale = from@parscale,
             random = from@random)
         # Note: groups slot is not copied since iwp doesn't have hierarchical structure
       })

#' Update model functions to return class objects
#'
#' These functions wrap the original model functions to return proper class objects
#' instead of lists.
#'
#' @name model_constructors
NULL
