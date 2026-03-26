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
           p = "numeric",
           ref_value = "numeric",
           knots = "ANY",
           range = "ANY",
           init = "numeric",
           lower = "numeric",
           upper = "numeric",
           parscale = "numeric"
         ),
         prototype = list(
           var = character(0),
           p = numeric(0),
           ref_value = numeric(0),
           knots = NULL,
           range = NULL,
           init = numeric(0),
           lower = numeric(0),
           upper = numeric(0),
           parscale = numeric(0)
         )
)

#' @rdname model_classes
#' @export
setClass("hiwp",
         representation = representation(
           var = "character",
           p = "numeric",
           ref_value = "numeric",
           knots = "ANY",
           range = "ANY",
           group_var = "character",
           init = "numeric",
           lower = "numeric",
           upper = "numeric",
           parscale = "numeric"
         ),
         prototype = list(
           var = character(0),
           p = numeric(0),
           ref_value = numeric(0),
           knots = NULL,
           range = NULL,
           group_var = character(0),
           init = numeric(0),
           lower = numeric(0),
           upper = numeric(0),
           parscale = numeric(0)
         )
)

#' @rdname model_classes
#' @export
setClass("fpoly",
         representation = representation(
           var = "character",
           p = "numeric",
           ref_value = "numeric",
           init = "numeric",
           lower = "numeric",
           upper = "numeric",
           parscale = "numeric"
         ),
         prototype = list(
           var = character(0),
           p = numeric(0),
           ref_value = numeric(0),
           init = numeric(0),
           lower = numeric(0),
           upper = numeric(0),
           parscale = numeric(0)
         )
)

#' @rdname model_classes
#' @export
setClass("rpoly",
         representation = representation(
           var = "character",
           p = "numeric",
           ref_value = "numeric",
           sd = "numeric"
         ),
         prototype = list(
           var = character(0),
           p = numeric(0),
           ref_value = numeric(0),
           sd = numeric(0)
         )
)

#' @rdname model_classes
#' @export
setClass("hrpoly",
         representation = representation(
           var = "character",
           p = "numeric",
           ref_value = "numeric",
           group_var = "character",
           init = "numeric",
           lower = "numeric",
           upper = "numeric",
           parscale = "numeric"
         ),
         prototype = list(
           var = character(0),
           p = numeric(0),
           ref_value = numeric(0),
           group_var = character(0),
           init = numeric(0),
           lower = numeric(0),
           upper = numeric(0),
           parscale = numeric(0)
         )
)

#' @rdname model_classes
#' @export
setClass("iid",
         representation = representation(
           var = "character",
           init = "numeric",
           lower = "numeric",
           upper = "numeric",
           parscale = "numeric"
         ),
         prototype = list(
           var = character(0),
           init = numeric(0),
           lower = numeric(0),
           upper = numeric(0),
           parscale = numeric(0)
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
             p = from@p,
             ref_value = from@ref_value,
             knots = from@knots,
             range = from@range,
             init = from@init,
             lower = from@lower,
             upper = from@upper,
             parscale = from@parscale)
       })

#' Update model functions to return class objects
#'
#' These functions wrap the original model functions to return proper class objects
#' instead of lists.
#'
#' @name model_constructors
NULL
