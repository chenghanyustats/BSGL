\name{BSGL-supporting-utilities}
\alias{BSGL-supporting-utilities}
\alias{beta_mse}
\alias{betas_bsgl}
\alias{betas_gs}
\alias{calc_cluster_beta}
\alias{calc_f1_beta}
\alias{calc_mse}
\alias{calc_scp}
\alias{calc_true_beta}
\alias{combine_chains}
\alias{compute_rhat}
\alias{compute_rhat_single}
\alias{cv_bsgl}
\alias{eval_bsgl}
\alias{eval_gs}
\alias{fit_bsgl_multichain}
\alias{fit_final_bsgl}
\alias{gen_data}
\alias{gen_rep}
\alias{get_gam_betas}
\alias{get_mgwr_betas}
\alias{load_data}
\alias{load_data_cluster}
\alias{make_basis}
\alias{plot_beta}
\alias{plot_coefficient_maps}
\alias{plot_observed_evi}
\alias{plot_observed_vs_predicted}
\alias{prepare_geographic_features}
\alias{prop_maps}
\alias{save_incl0_maps}
\alias{sum_bsgl}
\alias{sum_gs}
\alias{tensor.prod.model.matrix}
\title{Supporting Utilities for BSGL Analyses}
\description{
Supporting functions used by the BSGL simulation, diagnostics,
evaluation, visualization, data-preparation, and multichain workflows.
The primary model-fitting and prediction functions are documented on
their own help pages.
}
\details{
These utilities are retained in the exported package interface to support
the reproducibility workflows accompanying the BSGL manuscript.

Coefficient-surface and basis utilities include \code{betas_bsgl},
\code{betas_gs}, \code{make_basis}, \code{calc_true_beta},
\code{calc_cluster_beta}, and \code{tensor.prod.model.matrix}.

Evaluation and diagnostic utilities include \code{beta_mse},
\code{calc_mse}, \code{calc_f1_beta}, \code{calc_scp},
\code{compute_rhat}, \code{compute_rhat_single}, \code{cv_bsgl},
\code{eval_bsgl}, \code{eval_gs}, \code{sum_bsgl}, and \code{sum_gs}.

Simulation, data-loading, and multichain utilities include
\code{gen_data}, \code{gen_rep}, \code{load_data},
\code{load_data_cluster}, \code{combine_chains},
\code{fit_bsgl_multichain}, and \code{fit_final_bsgl}.

Comparison-model and visualization utilities include \code{get_gam_betas},
\code{get_mgwr_betas}, \code{plot_beta}, \code{plot_coefficient_maps},
\code{plot_observed_evi}, \code{plot_observed_vs_predicted},
\code{prepare_geographic_features}, \code{prop_maps}, and
\code{save_incl0_maps}.

For the main reviewer workflow, users will typically call the separately
documented model functions such as \code{fit_bsgl}, \code{pred_bsgl},
\code{fit_gs}, \code{pred_gs}, \code{fit_gam}, \code{pred_gam},
\code{fit_mgwr}, and \code{pred_mgwr} through the supplied demo scripts.
}
\value{
Return values depend on the utility. Functions return objects such as
numeric summaries, coefficient-surface matrices, fitted-model objects,
diagnostic summaries, plotting objects, or prepared data objects.
}
\seealso{
\code{fit_bsgl}, \code{pred_bsgl}, \code{fit_gs}, \code{pred_gs},
\code{fit_gam}, \code{pred_gam}, \code{fit_mgwr}, \code{pred_mgwr}
}
\keyword{utilities}
