

###  Association analysis between somatic units and OS


#' This function performs survival analysis for mutation counts mapped to PIU, LU and NCU.
#' 
#' @param piu_filename Filename for the PIU mapping results.
#' @param lu_filename Filename for the LU mapping results. 
#' @param ncu_filename Filename for the NCU mapping results. 
#' @param clinical_df Clinical information data frame. 
#' @param gender_as_covariate  Boolean variable indicating gender should be taken as a covariable. 
#' @param race_group_min The minimum number of patients with the same race enough for a race group.
#' @param min_surv_days The minimum number of survival time (in days) a patient has to be included in analysis.
#' @param min_surv_people The minimum number of patients survived (or censored) in a cohort that a survival analysis should be conducted.
#' @param patient_sum_min The minimum number of patients having mutation mapping to the unit for the unit to be analysed.
#' @param mutation_type A tag indicating the type of mutation; somatic or germline.
#' @param output_dir The directory you would like to have your output files in.
#' @import dplyr magrittr data.table qvalue survival
#' @export
#' @details 
#' @examples 
#' univariate_cox_model(piu_filename = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/example/piu_mapping_count.tsv",
#'                      lu_filename = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/example/lu_summarising_count.tsv",
#'                      ncu_filename = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/example/ncu_summarising_count.tsv",
#'                      clinical_df = sel_example_cdr,
#'                      gender_as_covariate = T,
#'                      race_group_min = 6,
#'                      min_surv_days = 90,
#'                      min_surv_people = 5,
#'                      patient_sum_min = 3,
#'                      mutation_type = "somatic",
#'                      output_dir = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/example/")
#'                       


univariate_cox_model = function(piu_filename,
                                lu_filename,
                                ncu_filename,
                                clinical_df,
                                gender_as_covariate = T,
                                race_group_min = 6,
                                min_surv_days = 90,
                                min_surv_people = 5,
                                patient_sum_min = 3,
                                mutation_type = "somatic",
                                output_dir)
  
{
  
  
  
  #Build cox proportional hazard model for PIUs
  
  univariate_cox_model_for_piu (piu_filename = piu_filename,
                                cdr_clinical = clinical_df,
                                gender_as_covariate = gender_as_covariate,
                                race_group_min = race_group_min,
                                min_surv_days = min_surv_days,
                                min_surv_people = min_surv_people,
                                patient_sum_min = patient_sum_min,
                                mutation_type = mutation_type,
                                output_dir = output_dir)
  cat("Models built for PIU.","\n")
  
  
  #Build cox proportional hazard model for LUs
  
  univariate_cox_model_for_lu (lu_filename = lu_filename,
                                 cdr_clinical = clinical_df,
                                 gender_as_covariate = gender_as_covariate,
                                 race_group_min = race_group_min,
                                 min_surv_days = min_surv_days,
                                 min_surv_people = min_surv_people,
                                 patient_sum_min = patient_sum_min,
                                 mutation_type = mutation_type,
                                 output_dir = output_dir)
  cat("Models built for LU.","\n")
  
  #Build cox proportional hazard model for NCUs
  
  univariate_cox_model_for_ncu (ncu_filename = ncu_filename,
                                cdr_clinical = clinical_df,
                                gender_as_covariate = gender_as_covariate,
                                race_group_min = race_group_min,
                                min_surv_days = min_surv_days,
                                min_surv_people = min_surv_people,
                                patient_sum_min = patient_sum_min,
                                mutation_type = mutation_type,
                                output_dir = output_dir)
  
  cat("Models built for NCU.","\n")
  
  
  
}
  
  


  