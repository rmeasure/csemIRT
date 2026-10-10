#' Simulated Attitude Scale Data
#'
#' Simulated responses to a 10-item, 5-point attitude scale for 500
#' respondents, generated from the Graded Response Model with slopes between
#' 1.3 and 2.0 and evenly spaced thresholds centred near zero.
#'
#' @format A data frame with 500 rows and 11 columns: item responses
#'   \code{A01} to \code{A10} (1-5) and \code{group}, a randomly assigned
#'   grouping variable ("Group_1", "Group_2") for subgroup analyses.
#' @source Simulated; see \code{data-raw/make_data.R} in the package sources.
#' @examples
#' data(attitude_data)
#' head(attitude_data)
"attitude_data"

#' Simulated Character Assessment Data
#'
#' Simulated responses to a 12-item, 5-point character assessment with
#' items labelled by content domain (cognitive, affective, behavioural) for
#' 400 respondents. Responses were generated from a unidimensional Graded
#' Response Model; the domain labels are descriptive only.
#'
#' @format A data frame with 400 rows and 13 columns: cognitive-domain
#'   items \code{C01}-\code{C04}, affective-domain items \code{A01}-\code{A04},
#'   behavioural-domain items \code{B01}-\code{B04} (all 1-5), and a randomly
#'   assigned \code{group} variable.
#' @source Simulated; see \code{data-raw/make_data.R}.
#' @examples
#' data(character_data)
#' head(character_data)
"character_data"

#' Simulated Well-Being Scale Data (GPCM)
#'
#' Simulated responses to a 12-item, 4-point positively worded well-being
#' scale for 600 respondents, generated from the Generalized Partial Credit
#' Model with strongly heterogeneous slopes (0.6 to 2.6) and step parameters
#' concentrated below the mean, producing a ceiling effect. The data are
#' designed to show how the choice among the GRM, GPCM and PCM changes the
#' estimated CSEM profile.
#'
#' @format A data frame with 600 rows and 13 columns: item responses
#'   \code{W01} to \code{W12} (1-4) and a randomly assigned \code{group}
#'   variable.
#' @source Simulated; see \code{data-raw/make_data.R}.
#' @examples
#' data(wellbeing_data)
#' csem_compare(wellbeing_data[, 1:12])
"wellbeing_data"
