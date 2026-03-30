# Helper function to convert hiwp to iwp
convert_hiwp_to_iwp <- function(from) {
  # Create a new iwp object with the same basic properties
  methods::new("iwp",
    term = from@term,
    formula = from@formula,
    knots = from@knots,
    ref_value = from@ref_value,
    p.order = from@p.order,
    by = character(0),  # hiwp has hierarchical structure, iwp doesn't
    init = from@init,
    lower = from@lower,
    upper = from@upper,
    parscale = from@parscale
  )
}

# Register the coercion method properly
methods::setAs("hiwp", "iwp", convert_hiwp_to_iwp)
methods::setAs("hrpoly", "rpoly", function(from) {
  rpoly(
    x = from@term,
    p = from@p.order,
    ref_value = from@ref_value,
    sd = 1
  )
})
