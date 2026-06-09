#' Evaluate TP and TN for Interaction Identification (CS-level evaluation)
#'
#' @param true_int_nam Character vector of true interactions, e.g. c("rs12345*Age", "rs45678*rs45678").
#' @param main_index Data frame from Identifying_MainEffect(): columns are Index, Variable (original names), and CS.
#' @param int_index Data frame from Identifying_IntEffect(): columns are Index, Variable (format "Main_CS1*Main_CS2" or "PM25*Main_CS1"), and CS.
#'
#' @return A list with two elements:
#'   - tp: 1 if all true interactions are identified; otherwise, 0.
#'   - tn: 1 if no false-positive interactions are identified; otherwise, 0.
#'
#' @export
int_tptn <- function(true_int_nam, main_index, int_index) {

# Step 1: Create mappings from CS to variable names, and from variable names to CS.
cs_to_var <- split(main_index$Variable, main_index$CS)
var_to_cs <- setNames(main_index$CS, main_index$Variable)

# Step 2: Function to expand interactions from CS-level into individual variable pairs.
expand_cs_pair <- function(part1, part2) {
  vars1 <- if(grepl("^Main_CS", part1)) cs_to_var[[part1]] else part1
  vars2 <- if(grepl("^Main_CS", part2)) cs_to_var[[part2]] else part2

  if(is.null(vars1) || is.null(vars2)) return(character(0))

  unique(c(
    paste0(rep(vars1, each=length(vars2)), "*", rep(vars2, times=length(vars1))),
    paste0(rep(vars2, each=length(vars1)), "*", rep(vars1, times=length(vars2)))
  ))
}

# Step 3: Identify CS pairs supported by true interactions.
supported_cs_pairs <- c()
for (pair in true_int_nam) {
  vars <- unlist(strsplit(pair, "\\*"))

  # Determine corresponding CS; if environment variable (e.g., PM25), use directly.
  cs1 <- if(vars[1] %in% names(var_to_cs)) var_to_cs[vars[1]] else vars[1]
  cs2 <- if(vars[2] %in% names(var_to_cs)) var_to_cs[vars[2]] else vars[2]

  # Add both directions for symmetry.
  supported_cs_pairs <- c(supported_cs_pairs,
                          paste(cs1, cs2, sep="*"),
                          paste(cs2, cs1, sep="*"))
}

# Step 4: Identify unsupported CS pairs (possible sources of false positives).
unsupported_cs_pairs <- setdiff(int_index$Variable, supported_cs_pairs)

# Step 5: Expand unsupported CS pairs into individual interactions (potential false positives).
risky_interactions <- unlist(lapply(strsplit(unsupported_cs_pairs, "\\*"),
                                    function(cs) expand_cs_pair(cs[1], cs[2])))

# Step 6: Expand all detected interactions for TP calculation.
detected <- unlist(lapply(strsplit(int_index$Variable, "\\*"),
                          function(cs) expand_cs_pair(cs[1], cs[2])))

# Step 7: Evaluate TP (true positive) and TN (true negative).
# TP: all true interactions must be detected.
tp <- as.integer(all(true_int_nam %in% detected))

# TN: no risky (unsupported) interactions should be detected.
tn <- as.integer(length(setdiff(risky_interactions, true_int_nam)) == 0)

return(list(tp=tp, tn=tn))
}
