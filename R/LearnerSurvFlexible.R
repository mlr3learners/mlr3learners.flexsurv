#' @title Survival Flexible Parametric Spline Learner
#'
#' @name mlr_learners_surv.flexible
#'
#' @description
#' A [mlr3proba::LearnerSurv] implementing flexible from package
#'   \CRANpkg{flexsurv}.
#' Calls [flexsurv::flexsurvspline()].
#'
#' @details
#' The `distr` prediction is estimated using the fitted custom distributions
#' from [flexsurv::flexsurvspline()] and the estimated coefficients however the prediction takes
#' place in this package and not in \CRANpkg{flexsurv} for a much faster and more efficient
#' implementation.
#'
#' As flexible spline models estimate the baseline hazard as the intercept, the linear predictor,
#' `lp`, can be calculated as in the classical setting. i.e. For fitted coefficients,
#' \eqn{\beta = (\beta_0,...,\beta_P)}{\beta = (\beta0,...,\betaP)},
#' and covariates \eqn{X^T = (X_0,...,X_P)^T}{X^T = (X0,...,XP)^T}, where \eqn{X_0}{X0} is a column
#' of \eqn{1}s: \eqn{lp = \beta X}{lp = \betaX}.
#'
#' @section Custom mlr3 defaults:
#' - `k`:
#'   - Actual default: `0`
#'   - Adjusted default: `1`
#'   - Reason for change: The default value of `0` is equivalent to, and a much less efficient
#'   implementation of, [LearnerSurvParametric][mlr3learners.survival::LearnerSurvParametric].
#'
#' @templateVar id surv.flexible
#' @template section_dictionary_learner
#'
#' @references
#' Royston P, Parmar MKB (2002).
#' “Flexible parametric proportional-hazards and proportional-odds models for censored survival
#' data, with application to prognostic modelling and estimation of treatment effects.”
#' Statistics in Medicine, 21(15), 2175–2197.
#' doi: 10.1002/sim.1203.
#'
#' @template seealso_learner
#' @template example
#' @export
LearnerSurvFlexible = R6Class("LearnerSurvFlexible",
  inherit = LearnerSurv,

  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      ps = ParamSet$new(
        params = list(
          ParamUty$new(id = "bhazard", tags = "train"),
          ParamInt$new(id = "k", default = 0L, lower = 0L, tags = "train"),
          ParamUty$new(id = "knots", tags = "train"),
          ParamUty$new(id = "bknots", tags = "train"),
          ParamFct$new(
            id = "scale", default = "hazard",
            levels = c("hazard", "odds", "normal"), tags = "train"),
          ParamFct$new(
            id = "timescale", default = "log",
            levels = c("log", "identity"), tags = "train"),
          ParamUty$new(id = "inits", tags = "train"),
          ParamUty$new(id = "fixedpars", tags = "train"),
          ParamDbl$new(id = "cl", default = 0.95, lower = 0, upper = 1, tags = "train"),
          ParamInt$new(id = "maxiter", default = 30L, tags = c("train", "control")),
          ParamDbl$new(id = "rel.tolerance", default = 1e-09, tags = c("train", "control")),
          ParamDbl$new(id = "toler.chol", default = 1e-10, tags = c("train", "control")),
          ParamInt$new(id = "debug", default = 0, lower = 0, upper = 1,
                       tags = c("train", "control")),
          ParamInt$new(id = "outer.max", default = 10L, tags = c("train", "control"))
      ))

      # value of k is changed as the default is equivalent (and a much more inefficient)
      # implementation of `surv.parametric`
      ps$values = list(k = 1)

      super$initialize(
        id = "surv.flexible",
        packages = "flexsurv",
        feature_types = c("logical", "integer", "factor", "numeric"),
        predict_types = c("distr", "crank", "lp"),
        param_set = ps,
        properties = "weights",
        man = "mlr3learners.flexsurv::mlr_learners_surv.flexible"
      )
    }
  ),

  private = list(
    .train = function(task) {
      pars_ctrl = self$param_set$get_values(tags = "control")
      pars_train = self$param_set$get_values(tags = "train")
      pars_train = pars_train[pars_train %nin% pars_ctrl]
      pars_train$sr.control = mlr3misc::invoke(survival::survreg.control, .args = pars_ctrl) #nolint

      if ("weights" %in% task$properties) {
        pars_train$weights = task$weights$weight
      }

      mlr3misc::invoke(flexsurv::flexsurvspline,
        formula = task$formula(task$feature_names),
        data = task$data(), .args = pars_train)
    },

    .predict = function(task) {

      # As we are using a custom predict method the missing assertions are performed here manually
      # (as opposed to the automatic assertions that take place after prediction)
      if (any(is.na(data.frame(task$data(cols = task$feature_names))))) {
        mlr3misc::stopf(
          "Learner %s on task %s failed to predict: Missing values in new data
                     (line(s) %s)\n",
          self$id, task$id,
          paste0(which(is.na(data.frame(task$data(cols = task$feature_names)))),
            collapse = ", "))
      }

      pred = mlr3misc::invoke(predict_flexsurvreg, self$model, task)

      # crank is defined as the mean of the survival distribution
      PredictionSurv$new(task = task, distr = pred$distr, lp = pred$lp, crank = pred$lp)
    }
  )
)
