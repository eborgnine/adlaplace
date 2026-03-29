# Helper function to convert hiwp to iwp
convert_hiwp_to_iwp <- function(from) {
  # Create a new iwp object with the same basic properties
  new("iwp",
    term = from@term,
    formula = term@formula,
    knots = from@knots,
    ref_value = from@ref_value,
    p.order = from@p.order,
    by = character(0),  # hiwp has hierarchical structure, iwp doesn't
    init = from@init,
    lower = from@lower,
    upper = from@upper,
    parscale = from@parscale,
    type = factor("random", levels = .type_factor_levels)
  )
}

# Register the coercion method properly
setAs("hiwp", "iwp", convert_hiwp_to_iwp)