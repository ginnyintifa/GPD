


####


nonzero = function(x) sum(x != 0)

cdr_tidy_up_for_model = function(interest_variable_info, unite_data, race_group_min,
                                 output_dir,
                                 output_name)
{


  patient_count = unite_data %>%
    dplyr::select(Tumor_Sample_Barcode, age,gender,one_of(interest_variable_info))

  unite_data$race[grep("\\[", unite_data$race)] = "OTHER"

  os_data = unite_data %>%
    dplyr::select(Tumor_Sample_Barcode, race, OS, OS.time) %>%
    na.omit()
  os_race_freq = as.data.frame(table(os_data$race))
  os_race_minor = os_race_freq %>%
    dplyr::filter(Freq<race_group_min)
  os_data$race[which(os_data$race %in% os_race_minor$Var1)] = "OTHER"
  colnames(os_data)[which(colnames(os_data)=="race")] = "os_race"


  patient_count_survival = patient_count %>%
    dplyr::left_join(os_data, by = "Tumor_Sample_Barcode") %>%
    dplyr::select(Tumor_Sample_Barcode, age, gender,
                  os_race, OS, OS.time,
                  one_of(interest_variable_info))





  write.table(patient_count_survival, paste0(output_dir, output_name),
              quote = F, row.names = F, sep = "\t")

  return(patient_count_survival)

}

 

piu_counts_cdr_clinical_unite = function(piu_count_filename,
                                         cdr_clinical,
                                         patient_sum_min,
                                         output_dir,
                                         output_name)

{
  
  # 
  # piu_count_filename = piu_filename
  # cdr_clinical = this_cdr
  # patient_sum_min = patient_sum_min
  # output_dir = output_dir
  # output_name = paste0(mutation_type,"_piu_cdr_clinical_unite.tsv")

  piu_count_df = fread(piu_count_filename,
                       stringsAsFactors = F)

  col_seq = seq(1:ncol(piu_count_df))
  which_count = col_seq[-c(1:8,ncol(piu_count_df))]
  # which_count = grep("TCGA", colnames(piu_count_df))
  non_zero = apply(piu_count_df[,..which_count],1,nonzero)
  piu_count_fsel = piu_count_df[which(non_zero>=patient_sum_min),]


  piu_count_sel = piu_count_fsel%>%
    dplyr::mutate(piu_info = paste(uniprot_accession, start_position,end_position,unit_label, unit_name, gene_name, gene_id, sep = "_"))%>%
    dplyr::mutate(gene_info = paste(gene_name, gene_id,sep = "_"))%>%
    dplyr::select(uniprot_accession, start_position, end_position, center_position, unit_name, gene_name, gene_id, unit_label,gene_info, piu_info, row_sum,
                  everything())


  if(nrow(piu_count_sel)>0)
  {

    col_seq = seq(1:ncol(piu_count_sel))
    which_barcode = col_seq[-c(1:11)]

    #which_barcode = grep("TCGA", colnames(piu_count_sel))

    piu_matrix = t(as.matrix(piu_count_sel[,which_barcode]))

    piu_count = data.frame(Tumor_Sample_Barcode = colnames(piu_count_sel)[which_barcode],
                           piu_matrix,
                           stringsAsFactors = F)

    colnames(piu_count) = c("Tumor_Sample_Barcode", piu_count_sel$piu_info)
    rownames(piu_count) = NULL

    piu_gene_df = data.frame(piu_info = piu_count_sel$piu_info,
                             gene_info = piu_count_sel$gene_info,
                             stringsAsFactors = F)


    piu_clinical_unite_data = piu_count %>%
      dplyr::left_join(cdr_clinical, by = "Tumor_Sample_Barcode") %>%
      dplyr::select(Tumor_Sample_Barcode,age, gender, race, OS, OS.time,everything())


    return_list = vector(mode = "list", length = 4)

    return_list[[1]] = piu_clinical_unite_data
    return_list[[2]] = piu_gene_df
    return_list[[3]] = piu_count_sel$piu_info
    return_list[[4]] = piu_count_sel$gene_id




    return(return_list)

  }else{
    cat("No item on this type of PIU.","\n")
  }


}



gene_counts_cdr_clinical_unite = function(gene_count_filename,
                                          cdr_clinical,
                                          patient_sum_min,
                                          output_dir,
                                          output_name)

{

  gene_count_df = fread(gene_count_filename,
                        stringsAsFactors = F)


  if("gene_info" %in% colnames(gene_count_df))
  {
    gene_count_df_info = as_tibble(gene_count_df)%>%
      dplyr::select(gene_id, gene_name, gene_info, everything())

  }else{
    gene_count_df_info = gene_count_df %>%
      dplyr::mutate(gene_info = paste(gene_name, gene_id, sep = "_"))%>%
      dplyr::select(gene_id, gene_name, gene_info, everything())

  }


  col_seq = seq(1:ncol(gene_count_df_info))
  which_count = col_seq[-c(1:3,ncol(gene_count_df_info))]

  # which_count = grep("TCGA", colnames(gene_count_df_info))
  non_zero = apply(gene_count_df_info[,which_count],1,nonzero)
  gene_count_sel = gene_count_df_info[which(non_zero>=patient_sum_min),]



  col_seq = seq(1:ncol(gene_count_sel))
  which_barcode = col_seq[-c(1:3,ncol(gene_count_sel))]
  #which_barcode = grep("TCGA", colnames(gene_count_sel))

  gene_matrix = t(as.matrix(gene_count_sel[,which_barcode]))

  gene_count = data.frame(Tumor_Sample_Barcode = colnames(gene_count_sel)[which_barcode],
                          gene_matrix,
                          stringsAsFactors = F)

  colnames(gene_count) = c("Tumor_Sample_Barcode", gene_count_sel$gene_info)
  rownames(gene_count) = NULL


  gene_clinical_unite_data = gene_count %>%
    dplyr::left_join(cdr_clinical, by = "Tumor_Sample_Barcode") %>%
    dplyr::select(Tumor_Sample_Barcode,age, gender, race, OS, OS.time,everything())


  return_list = vector(mode = "list", length = 2)

  return_list[[1]] = gene_clinical_unite_data
  return_list[[2]] = gene_count_sel$gene_info



  return(return_list)

}






filter_survival_data = function (surv_info_data,
                                 surv_status = quo(OS),
                                 surv_time = quo(OS.time),
                                 surv_race = quo(os_race),
                                 min_surv_days = 90)
  
  
{
  # 
  # surv_info_data = survival_info[survival_train, ]
  # surv_status = quo(OS)
  # surv_time = quo(OS.time)
  # surv_race = quo(os_race)
  # min_surv_days = 90
  # 
  unit_names = grep("ENSG", colnames(surv_info_data), value = T)
  
  
  filter_surv_data = surv_info_data %>%
    dplyr::select(barcode, age,gender, !!surv_race, !!surv_status, !!surv_time, one_of(unit_names)) %>%
    na.omit()%>%
    # dplyr::arrange(!!surv_status, !!surv_time) %>%
    dplyr::filter(!!surv_time >= min_surv_days)%>%
    na.omit()
  
  colnames(filter_surv_data) = c("Tumor_Sample_Barcode", "age", "gender", "race", "survival_status","survival_time", unit_names)
  
  
  return(filter_surv_data)
  
  
}








univariate_survival_significance_adjust_totalMutStage = function(filter_surv_data,
                                                                  output_dir,
                                                                  surv_name)
{

  age_data = as.numeric(filter_surv_data$age)
  stage_data = as.numeric(filter_surv_data$code_stage)
  totalMut_data = as.numeric(filter_surv_data$total_mutation)
  gender_data = relevel(as.factor(filter_surv_data$gender), ref = unique(filter_surv_data$gender)[1])
  race_data = relevel(as.factor(filter_surv_data$race), ref = unique(filter_surv_data$race)[1])
  surv_object = Surv(time = filter_surv_data$survival_time, event = filter_surv_data$survival_status)


  unit_names = grep("ENSG", colnames(filter_surv_data), value = T)

  surv_result_df = rbindlist(lapply(1:length(unit_names), function(x)

    #for(x in 1:length(unit_names))
  {
    #cat(x, "\n")

    #x = 10

    this_count = filter_surv_data[[unit_names[x]]]

    if(length(which(this_count>0))>=2)   #### do not build model if the major effect is to sparse
    {

      this_mutate_patients = sum(this_count!= 0)


      this_count_coeff = NA
      this_count_exp_coeff = NA
      this_count_pval = NA
      this_assum_test = NA


      this_total_coeff = NA
      this_total_exp_coeff = NA
      this_total_pval = NA
      this_total_test = NA


      if(length(unique(gender_data))>1)
      {
        if(length(unique(race_data))>1)
        {
          this_model = coxph(surv_object ~  totalMut_data + age_data  + gender_data + race_data + stage_data +
                               this_count)
        }else{
          this_model = coxph(surv_object ~ totalMut_data + age_data  + gender_data + stage_data +
                               this_count)
        }

      }else{
        if(length(unique(race_data))>1)
        {
          this_model = coxph(surv_object ~  totalMut_data + age_data  + race_data + stage_data +
                               this_count)
        }else{
          this_model = coxph(surv_object ~  totalMut_data + age_data + stage_data +
                               this_count)
        }

      }


      this_test = cox.zph(this_model)
      this_table = this_test$table



      #### after lunch


      this = summary(this_model)
      this_coef = this$coefficients
      r_this =  length(unique(gender_data)) +length(unique(race_data)) +1+1
      r_stage =  length(unique(gender_data)) +length(unique(race_data)) +1

      this_count_coeff = this_coef[r_this,1]
      this_count_exp_coeff = this_coef[r_this,2]
      this_count_pval = this_coef[r_this,5]
      this_assum_test = this_table[r_this,3]


      this_total_coeff = this_coef[1,1]
      this_total_exp_coeff = this_coef[1,2]
      this_total_pval = this_coef[1,5]
      this_total_test = this_table[1,3]




      one_surv_result = data.frame(count_info = unit_names[x],
                                   mutation_patients = this_mutate_patients,
                                   total_patients = nrow(filter_surv_data),
                                   count_coeff = this_count_coeff,
                                   count_exp_coeff = this_count_exp_coeff,
                                   count_pval = this_count_pval,
                                   assum_test = this_assum_test,
                                   count_qval = NA,
                                   total_coeff = this_total_coeff,
                                   total_exp_coeff = this_total_exp_coeff,
                                   total_pval = this_total_pval,
                                   total_test = this_total_test,
                                   stage_coeff = this_coef[r_stage, 1],
                                   stage_exp_coeff = this_coef[r_stage,2],
                                   stage_count_pval = this_coef[r_stage, 5],
                                   stage_test = this_table[r_stage, 3],


                                   # total_qval = NA,
                                   stringsAsFactors = F)

      # if(x%%2000 == 0)
      # cat(x,"/", length(unit_names),  "\n")
      #
      return(one_surv_result)

    }




  }))



  if(length(unit_names)>5)
  {

    if(length(surv_result_df$count_pval)>300 & min(surv_result_df$count_pval,na.rm = T)<0.05 & max(surv_result_df$count_pval,na.rm = T)>0.95)
    {
      qval = qvalue(surv_result_df$count_pval)
    }else{
      qval = qvalue(surv_result_df$count_pval, pi0 = 1)

    }
    surv_result_df$count_qval = qval$qvalues



  }else{
    cat("Small unit size, no q-val estimation.", "\n")

  }



  output_name = paste0(surv_name,".tsv")

  return(surv_result_df)

}


##### running the following function requres two addition file describing the stage information and total mutation information.
####  Their names are barcode_stage_filename and somatic_mc3_filename
#### Format of barcode_stage_filename:
# Tumor_Sample_Barcode code_stage
# example-1          1
# example-2          2
# example-3          3
# example-4          4
# ...
### Format of somatic_mc3_filename:
# Hugo_Symbol     Gene    Chromosome      Start_Position  End_Position    Variant_Classification  Variant_Type    HGVSc   HGVSp   agg_sample_id   mut_freq
# CTNNB1  ENSG00000168036 3       41266136        41266136        Missense_Mutation       SNP     c.133T>C        p.Ser45Pro      example-1_example-2_example-3        3
# AC005013.5      ENSG00000228421 7       28996482        28996482        5'Flank DEL     .       .       example-2_example3       2
# DLGAP3  ENSG00000116544 1       35365799        35365799        Missense_Mutation       SNP     c.1183G>A       p.Gly395Ser     example-1_example2       2
# OR4K2   ENSG00000165762 14      20345046        20345046        Missense_Mutation       SNP     c.620C>T        p.Ala207Val     example-2_example4       2
# ...




univariate_cox_model_adjustTotalMutationStage = function(piu_filename,
                                lu_filename,
                                ncu_filename,
                                barcode_stage_filename,
                                somatic_mc3_filename,
                                clinical_df,
                                race_group_min = 6,
                                min_surv_days = 90,
                                min_surv_people = 5,
                                patient_sum_min = 3,
                                mutation_type = "somatic",
                                output_dir)

{

  get_barcode_stage = fread(barcode_stage_filename, stringsAsFactors = F)

  cancer_mc3 = fread(somatic_mc3_filename, stringsAsFactors = F)



  #### piu

  if(file.exists(piu_filename))
  {
    piu_unite = piu_counts_cdr_clinical_unite(
      piu_count_filename = piu_filename,
      cdr_clinical = cdr_clinical,
      patient_sum_min = patient_sum_min,
      output_dir = output_dir,
      output_name = paste0(mutation_type,"_piu_cdr_clinical_unite.tsv"))

    if(file.exists(paste0(output_dir,paste0(mutation_type,"_piu_cdr_clinical_unite.tsv"))))
    {
      surv_piu_info = cdr_tidy_up_for_model(
        interest_variable_info = piu_unite[[3]],
        unite_data = piu_unite[[1]],
        race_group_min = race_group_min,
        output_dir = output_dir,
        output_name = paste0(mutation_type,"_piu_survival_info.tsv"))

    }
  }

    ####bpiu

    if(file.exists(bpiu_filename))
    {

      ### an additional step for selecting valid genes


      bpiu_unite = gene_counts_cdr_clinical_unite(
        gene_count_filename = bpiu_filename,
        cdr_clinical = cdr_clinical,
        patient_sum_min = patient_sum_min,
        output_dir = output_dir,
        output_name = paste0(mutation_type,"_bpiu_cdr_clinical_unite.tsv"))



      if(file.exists(paste0(output_dir,paste0(mutation_type,"_bpiu_cdr_clinical_unite.tsv"))))
      {
        surv_bpiu_info = cdr_tidy_up_for_model(
          interest_variable_info = bpiu_unite[[2]],
          unite_data = bpiu_unite[[1]],
          race_group_min = race_group_min,
          output_dir = output_dir,
          output_name = paste0(mutation_type,"_bpiu_survival_info.tsv"))

      }
    }

  ####npc



  if(file.exists(npc_filename))
  {



    npc_unite = gene_counts_cdr_clinical_unite(
      gene_count_filename = npc_filename,
      cdr_clinical = cdr_clinical,
      patient_sum_min = patient_sum_min,
      output_dir = output_dir,
      output_name = paste0(mutation_type,"_npc_cdr_clinical_unite.tsv"))



    if(file.exists(paste0(output_dir,paste0(mutation_type,"_npc_cdr_clinical_unite.tsv"))))
    {
      surv_npc_info = cdr_tidy_up_for_model(
        interest_variable_info = npc_unite[[2]],
        unite_data = npc_unite[[1]],
        race_group_min = race_group_min,
        output_dir = output_dir,
        output_name = paste0(mutation_type,"_npc_survival_info.tsv"))


    }
  }



  surv_names = grep("ENSG", colnames(surv_piu_info), value = T, invert = T)

  piu_names = grep("ENSG", colnames(surv_piu_info), value = T)

  bpiu_names = grep("ENSG", colnames(surv_bpiu_info), value = T)

  npc_names = grep("ENSG", colnames(surv_npc_info), value = T)

  if(length(piu_names) >0)
  {
    colnames(surv_piu_info) = c(surv_names,paste0("_piu_", piu_names))

  }


  if(length(bpiu_names)>0)
  {
    colnames(surv_bpiu_info) = c(surv_names,paste0("_bpiu_", bpiu_names))

  }

  if(length(npc_names)>0)
  {

    colnames(surv_npc_info) = c(surv_names,paste0("_npc_", npc_names))

  }



  surv_info = surv_piu_info%>%
    dplyr::left_join(surv_bpiu_info, by = surv_names)%>%
    dplyr::left_join(surv_npc_info, by = surv_names)%>%
    dplyr::filter(barcode%in%cdr_no_hyper$bcr_patient_barcode)



  filter_surv_os_data = filter_survival_data(surv_info_data = surv_info,
                                             surv_status = quo(OS),
                                             surv_time = quo(OS.time),
                                             surv_race = quo(os_race),
                                             min_surv_days = 90)





  ns_mc3 = cancer_mc3%>%
    dplyr::filter(! Variant_Classification == "Silent")


  all_barcodes = unlist(lapply(1:nrow(ns_mc3), function(x) {


    this_row = ns_mc3$agg_sample_id[x]

    sep_id = unlist(strsplit(this_row, split = "_", fixed = T))



  }))



  barcode_df = as.data.frame(table(all_barcodes))%>%
    dplyr::arrange(desc(Freq))

  barcode_totalMut = data.frame(Tumor_Sample_Barcode = as.character(barcode_df$all_barcodes),
                                total_mutation = barcode_df$Freq,
                                stringsAsFactors = F)



  surv_with_totalMut = barcode_totalMut %>%
    dplyr::left_join(filter_surv_os_data, by = "Tumor_Sample_Barcode")%>%
    dplyr::left_join(get_barcode_stage, by = "Tumor_Sample_Barcode")%>%
    na.omit()

  test_units = univariate_survival_significance_adjust_totalMutStage(filter_surv_data = surv_with_totalMut,
                                                                output_dir = output_dir,
                                                                surv_name = "os_significance")



  write.table(test_units, paste0(output_dir,"all_units_results.tsv"),
              quote = F, row.names = F, sep = "\t")




}




