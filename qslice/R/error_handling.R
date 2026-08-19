## Error-handling framework suggested by ChatGPT Codex
## Implemented and modified by Matt Heiner

.format_state <- function(x, max_chars = 300L) {
  txt <- paste(capture.output(dput(x)), collapse = "")
  if (nchar(txt) > max_chars) {
    paste0(substr(txt, 1L, max_chars), "...")
  } else {
    txt
  }
}


.eval_log_target <- function(
  log_target,
  x,
  sampler,
  stage,
  evaluation,
  allow_minus_inf = TRUE
) {
  fail <- function(reason, hint, value = NULL, parent = NULL) {
    state_txt <- .format_state(x)

    value_txt <- if (is.null(value)) {
      "<log-target (or likelihood) evaluation not available>"
    } else {
      .format_state(value)
    }

    message <- paste0(
      sampler,
      ": log-target (or likelihood) evaluation failed.\n",
      "  Stage: ",
      stage,
      "\n",
      "  Evaluation: ",
      evaluation,
      "\n",
      "  Proposed/current state: ",
      state_txt,
      "\n",
      "  log-density evaluation: ",
      value_txt,
      "\n",
      "  Reason: ",
      reason,
      "\n",
      "  Likely explanation: ",
      hint
    )

    if (!is.null(parent)) {
      message <- paste0(
        message,
        "\n  Original log-target error: ",
        conditionMessage(parent)
      )
    }

    stop(errorCondition(
      message = message,
      class = c("qslice_target_error", "qslice_error"),
      sampler = sampler,
      stage = stage,
      evaluation = evaluation,
      state = x,
      target_value = value,
      reason = reason,
      parent = parent
    ))
  }

  # Diagnose a bad proposal before passing it to the target.
  if (!is.numeric(x) || length(x) < 1L || anyNA(x) || any(!is.finite(x))) {
    fail(
      reason = "The sampler generated an invalid state.",
      hint = paste0(
        "A proposal, quantile function, or numerical calculation produced ",
        "a nonnumeric, missing, or nonfinite value."
      )
    )
  }

  value <- tryCatch(
    log_target(x),
    error = function(e) {
      if (inherits(e, "qslice_error")) {
        stop(e)
      }
      fail(
        reason = "The log-target (or likelihood) function threw an R error.",
        hint = paste0(
          "Inspect the original error below. Common causes include an ",
          "invalid parameter, indexing error, or input outside the ",
          "function's intended domain."
        ),
        parent = e
      )
    }
  )

  if (!is.numeric(value) || length(value) != 1L) {
    fail(
      reason = sprintf(
        "The log-target (or likelihood) returned type `%s` with length %d.",
        typeof(value),
        length(value)
      ),
      hint = "A log-target must return exactly one numeric value.",
      value = value
    )
  }

  # Check NaN before NA because is.na(NaN) is also TRUE.
  if (is.nan(value)) {
    fail(
      reason = "The log-target (or likelihood) returned NaN.",
      hint = paste0(
        "NaN usually results from an undefined operation, such as 0/0, ",
        "Inf - Inf, log of a negative value, or an invalid distribution parameter."
      ),
      value = value
    )
  }

  if (is.na(value)) {
    fail(
      reason = "The log-target (or likelihood) returned NA.",
      hint = paste0(
        "An input or intermediate target calculation is missing, or a ",
        "conditional branch did not return a valid value."
      ),
      value = value
    )
  }

  if (value == Inf) {
    fail(
      reason = "The log-target (or likelihood) returned +Inf.",
      hint = paste0(
        "This commonly indicates numerical overflow, a singularity, or ",
        "an incorrectly specified log density."
      ),
      value = value
    )
  }

  if (value == -Inf && !allow_minus_inf) {
    fail(
      reason = "The current state has zero target (or likelihood) density.",
      hint = paste0(
        "The initial/current state is outside the target support. ",
        "Initialize the chain at a state with finite log density."
      ),
      value = value
    )
  }

  value
}


.eval_pseudo_logdens <- function(
  logdens,
  x,
  sampler,
  stage,
  evaluation,
  allow_minus_inf = FALSE
) {
  fail <- function(reason, hint, value = NULL, parent = NULL) {
    state_txt <- .format_state(x)

    value_txt <- if (is.null(value)) {
      "<log-density evaluation not available>"
    } else {
      .format_state(value)
    }

    message <- paste0(
      sampler,
      ": pseudo-target log-density evaluation failed.\n",
      "  Stage: ",
      stage,
      "\n",
      "  Evaluation: ",
      evaluation,
      "\n",
      "  Proposed/current state: ",
      state_txt,
      "\n",
      "  log-density evaluation: ",
      value_txt,
      "\n",
      "  Reason: ",
      reason,
      "\n",
      "  Likely explanation: ",
      hint
    )

    if (!is.null(parent)) {
      message <- paste0(
        message,
        "\n  Original pseudo-target log-density error: ",
        conditionMessage(parent)
      )
    }

    stop(errorCondition(
      message = message,
      class = c("qslice_pseudo_density_error", "qslice_error"),
      sampler = sampler,
      stage = stage,
      evaluation = evaluation,
      state = x,
      logdensity_value = value,
      reason = reason,
      parent = parent
    ))
  }

  # Diagnose a bad proposal before passing it.
  if (!is.numeric(x) || length(x) < 1L || anyNA(x) || any(!is.finite(x))) {
    fail(
      reason = "The sampler generated an invalid state.",
      hint = paste0(
        "A proposal, quantile function, or numerical calculation produced ",
        "a nonnumeric, missing, or nonfinite value."
      )
    )
  }

  value <- tryCatch(
    logdens(x),
    error = function(e) {
      if (inherits(e, "qslice_error")) {
        stop(e)
      }
      fail(
        reason = "The pseudo-target log-density function threw an R error.",
        hint = paste0(
          "Inspect the original error below. Common causes include an ",
          "invalid parameter, indexing error, or input outside the ",
          "function's intended domain."
        ),
        parent = e
      )
    }
  )

  if (!is.numeric(value) || length(value) != 1L) {
    fail(
      reason = sprintf(
        "The pseudo-target log-density returned type `%s` with length %d.",
        typeof(value),
        length(value)
      ),
      hint = "A pseudo-target log-density must return exactly one numeric value.",
      value = value
    )
  }

  # Check NaN before NA because is.na(NaN) is also TRUE.
  if (is.nan(value)) {
    fail(
      reason = "The pseudo-target log-density returned NaN.",
      hint = paste0(
        "NaN usually results from an undefined operation, such as 0/0, ",
        "Inf - Inf, log of a negative value, or an invalid distribution parameter."
      ),
      value = value
    )
  }

  if (is.na(value)) {
    fail(
      reason = "The pseudo-target log-density returned NA.",
      hint = paste0(
        "An input or intermediate target calculation is missing, or a ",
        "conditional branch did not return a valid value."
      ),
      value = value
    )
  }

  if (value == Inf) {
    fail(
      reason = "The pseudo-target log-density returned +Inf.",
      hint = paste0(
        "This commonly indicates numerical overflow, a singularity, or ",
        "an incorrectly specified log density."
      ),
      value = value
    )
  }

  if (value == -Inf && !allow_minus_inf) {
    fail(
      reason = "The current state has zero pseudo-target density.",
      hint = paste0(
        "The initial/current state is outside the pseudo-target support. ",
        "This may result from changing the pseudo-target to one that does not support the most recent chains state."
      ),
      value = value
    )
  }

  value
}

.eval_pseudo_quantile <- function(
  q,
  u,
  sampler,
  stage,
  attempt,
  expected_length,
  require_interior = TRUE
) {
  fail <- function(
    reason,
    hint,
    value = NULL,
    parent = NULL,
    warnings = character()
  ) {
    message <- paste0(
      sampler,
      ": pseudo-target quantile evaluation failed.\n",
      "  Stage: ",
      stage,
      "\n",
      "  Proposal attempt: ",
      attempt,
      "\n",
      "  Quantile probability u: ",
      .format_state(u),
      "\n",
      "  pseudo$q(u): ",
      if (is.null(value)) "<not available>" else .format_state(value),
      "\n",
      "  Reason: ",
      reason,
      "\n",
      "  Likely explanation: ",
      hint
    )

    if (length(warnings)) {
      message <- paste0(
        message,
        "\n  Warnings from pseudo$q(): ",
        paste(unique(warnings), collapse = "; ")
      )
    }

    if (!is.null(parent)) {
      message <- paste0(
        message,
        "\n  Original pseudo$q() error: ",
        conditionMessage(parent)
      )
    }

    stop(errorCondition(
      message = message,
      class = c("qslice_pseudo_quantile_error", "qslice_error"),
      sampler = sampler,
      stage = stage,
      attempt = attempt,
      probability = u,
      quantile_value = value,
      reason = reason,
      warnings = warnings,
      parent = parent
    ))
  }

  if (!is.function(q)) {
    fail(
      reason = "`pseudo$q` is not a function.",
      hint = paste0(
        "The pseudo-target must contain a quantile function named `q`, ",
        "for example `list(ld = ..., q = ...)`."
      )
    )
  }

  if (!is.numeric(u) || length(u) < 1L) {
    fail(
      reason = "`u` is not a nonempty numeric value.",
      hint = "The sampler generated an invalid quantile probability."
    )
  }

  if (anyNA(u) || any(!is.finite(u))) {
    fail(
      reason = "`u` contains NA, NaN, or an infinite value.",
      hint = paste0(
        "The shrinking interval or a preceding pseudo-target calculation ",
        "produced an invalid probability."
      )
    )
  }

  if (any(u < 0 | u > 1)) {
    fail(
      reason = "`u` lies outside [0, 1].",
      hint = paste0(
        "This indicates an error in the shrinking interval or in a ",
        "pseudo-target CDF calculation."
      )
    )
  }

  if (require_interior && any(u <= 0 | u >= 1)) {
    fail(
      reason = "`u` is not strictly inside (0, 1).",
      hint = paste0(
        "Quantile-slice proposals should use interior probabilities. ",
        "An exact endpoint can produce an infinite quantile for an ",
        "unbounded pseudo-target."
      )
    )
  }

  q_warnings <- character()

  value <- tryCatch(
    withCallingHandlers(
      q(u),
      warning = function(w) {
        q_warnings <<- c(q_warnings, conditionMessage(w))
      }
    ),
    error = function(e) {
      if (inherits(e, "qslice_error")) {
        stop(e)
      }
      fail(
        reason = "`pseudo$q(u)` threw an R error.",
        hint = paste0(
          "The quantile function may have invalid parameters, may not ",
          "accept this input shape, or may contain an implementation error."
        ),
        parent = e,
        warnings = q_warnings
      )
    }
  )

  if (!is.numeric(value)) {
    fail(
      reason = sprintf(
        "`pseudo$q(u)` returned an object of type `%s`.",
        typeof(value)
      ),
      hint = "The quantile function must return numeric proposal values.",
      value = value,
      warnings = q_warnings
    )
  }

  if (length(value) != expected_length) {
    fail(
      reason = sprintf(
        "`pseudo$q(u)` returned length %d; expected length %d.",
        length(value),
        expected_length
      ),
      hint = paste0(
        "Check whether the quantile function is properly vectorized and ",
        "whether it returns one proposal per requested component."
      ),
      value = value,
      warnings = q_warnings
    )
  }

  if (any(is.nan(value))) {
    fail(
      reason = "`pseudo$q(u)` returned NaN.",
      hint = paste0(
        "This usually indicates invalid distribution parameters or ",
        "numerical failure in an extreme tail."
      ),
      value = value,
      warnings = q_warnings
    )
  }

  if (any(is.na(value))) {
    fail(
      reason = "`pseudo$q(u)` returned NA.",
      hint = paste0(
        "The quantile function may have received invalid parameters or ",
        "propagated a missing intermediate value."
      ),
      value = value,
      warnings = q_warnings
    )
  }

  if (any(!is.finite(value))) {
    location <- if (any(u <= 0 | u >= 1)) {
      "at a boundary probability"
    } else {
      "at an interior probability"
    }

    fail(
      reason = paste0(
        "`pseudo$q(u)` returned an infinite proposal ",
        location,
        "."
      ),
      hint = if (any(u > 0 & u < 1)) {
        paste0(
          "An infinite quantile at an interior probability usually means ",
          "the pseudo-target parameters are invalid or its tail calculation ",
          "is numerically unstable."
        )
      } else {
        paste0(
          "Infinite endpoint quantiles are mathematically valid for ",
          "unbounded distributions, but they are not valid sampler states."
        )
      },
      value = value,
      warnings = q_warnings
    )
  }

  value
}

.eval_pseudo_cdf <- function(
  p,
  x,
  sampler,
  stage,
  attempt,
  expected_length,
  require_interior = FALSE
) {
  fail <- function(
    reason,
    hint,
    value = NULL,
    parent = NULL,
    warnings = character()
  ) {
    message <- paste0(
      sampler,
      ": pseudo-target CDF evaluation failed.\n",
      "  Stage: ",
      stage,
      "\n",
      "  Proposal attempt: ",
      attempt,
      "\n",
      "  State supplied to pseudo$p(): ",
      .format_state(x),
      "\n",
      "  pseudo$p(x): ",
      if (is.null(value)) "<not available>" else .format_state(value),
      "\n",
      "  Reason: ",
      reason,
      "\n",
      "  Likely explanation: ",
      hint
    )

    if (length(warnings)) {
      message <- paste0(
        message,
        "\n  Warnings from pseudo$p(): ",
        paste(unique(warnings), collapse = "; ")
      )
    }

    if (!is.null(parent)) {
      message <- paste0(
        message,
        "\n  Original pseudo$p() error: ",
        conditionMessage(parent)
      )
    }

    stop(errorCondition(
      message = message,
      class = c("qslice_pseudo_cdf_error", "qslice_error"),
      sampler = sampler,
      stage = stage,
      attempt = attempt,
      state = x,
      cdf_value = value,
      reason = reason,
      warnings = warnings,
      parent = parent
    ))
  }

  if (!is.function(p)) {
    fail(
      reason = "`pseudo$p` is not a function.",
      hint = paste0(
        "The pseudo-target must contain a CDF function named `p`, ",
        "for example `list(ld = ..., p = ..., q = ...)`."
      )
    )
  }

  if (!is.numeric(x) || is.complex(x) || length(x) < 1L) {
    fail(
      reason = "`x` is not a nonempty real numeric value.",
      hint = "The sampler passed an invalid state to the pseudo-target CDF."
    )
  }

  if (anyNA(x) || any(!is.finite(x))) {
    fail(
      reason = "`x` contains NA, NaN, or an infinite value.",
      hint = paste0(
        "A proposal, quantile calculation, or preceding sampler step ",
        "produced an invalid state."
      )
    )
  }

  p_warnings <- character()

  value <- tryCatch(
    withCallingHandlers(
      p(x),
      warning = function(w) {
        p_warnings <<- c(p_warnings, conditionMessage(w))
      }
    ),
    error = function(e) {
      if (inherits(e, "qslice_error")) {
        stop(e)
      }
      fail(
        reason = "`pseudo$p(x)` threw an R error.",
        hint = paste0(
          "The CDF may have invalid parameters, may not accept this ",
          "input shape, or may contain an implementation error."
        ),
        parent = e,
        warnings = p_warnings
      )
    }
  )

  if (!is.numeric(value) || is.complex(value)) {
    fail(
      reason = sprintf(
        "`pseudo$p(x)` returned an object of type `%s`.",
        typeof(value)
      ),
      hint = "The CDF must return real numeric probabilities.",
      value = value,
      warnings = p_warnings
    )
  }

  if (length(value) != expected_length) {
    fail(
      reason = sprintf(
        "`pseudo$p(x)` returned length %d; expected length %d.",
        length(value),
        expected_length
      ),
      hint = paste0(
        "Check whether the CDF is properly vectorized and whether it ",
        "returns one probability per supplied state."
      ),
      value = value,
      warnings = p_warnings
    )
  }

  # Check NaN before NA because is.na(NaN) is also TRUE.
  if (any(is.nan(value))) {
    fail(
      reason = "`pseudo$p(x)` returned NaN.",
      hint = paste0(
        "This usually indicates invalid distribution parameters, an ",
        "undefined calculation, or numerical failure in an extreme tail."
      ),
      value = value,
      warnings = p_warnings
    )
  }

  if (any(is.na(value))) {
    fail(
      reason = "`pseudo$p(x)` returned NA.",
      hint = paste0(
        "The CDF may have received invalid parameters or propagated ",
        "a missing intermediate value."
      ),
      value = value,
      warnings = p_warnings
    )
  }

  if (any(!is.finite(value))) {
    fail(
      reason = "`pseudo$p(x)` returned an infinite value.",
      hint = paste0(
        "A CDF must return a finite probability between zero and one. ",
        "Check the pseudo-target parameters and its CDF implementation."
      ),
      value = value,
      warnings = p_warnings
    )
  }

  if (any(value < 0 | value > 1)) {
    fail(
      reason = "`pseudo$p(x)` returned a value outside [0, 1].",
      hint = paste0(
        "The function supplied as `pseudo$p` is not producing a valid ",
        "CDF value, possibly because of invalid parameters or a coding error."
      ),
      value = value,
      warnings = p_warnings
    )
  }

  if (require_interior && any(value <= 0 | value >= 1)) {
    fail(
      reason = "`pseudo$p(x)` returned a boundary probability.",
      hint = paste0(
        "The sampler requires a transformed state strictly inside (0, 1). ",
        "The state may lie outside the pseudo-target support, or the CDF ",
        "may have numerically underflowed to 0 or rounded to 1 in an ",
        "extreme tail."
      ),
      value = value,
      warnings = p_warnings
    )
  }

  value
}
