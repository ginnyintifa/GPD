


univariate_cox_model_for_piu = function(piu_filename,
                                                cdr_clinical,
                                                gender_as_covariate = T,
                                                race_group_min = 6,
                                                min_surv_days = 90,
                                                min_surv_people = 5,
                                                patient_sum_min = 3,
                                                mutation_type = "somatic",
                                                output_dir)
{
  
   
  if(file.exists(piu_filename))
  {
    piu_unite = piu_counts_cdr_clinical_unite(
      piu_count_filename = piu_filename,
      cdr_clinical = cdr_clinical,
      patient_sum_min = patient_sum_min,
      output_dir = output_dir,
      output_name = paste0(mutation_type,"_piu_cdr_clinical_unite.tsv"))
    
    if(nrow(piu_unite[[1]])>0)
    {
      piu_info = cdr_tidy_up_for_model(
        interest_variable_info = piu_unite[[3]],
        unite_data = piu_unite[[1]],
        race_group_min = race_group_min,
        output_dir = output_dir,
        output_name = paste0(mutation_type,"_piu_survival_info.tsv"))
      
      if(gender_as_covariate == T)
      {
        survival_model(
          surv_info_data = piu_info,
          interest_variable_info = piu_unite[[3]],
          min_surv_time = min_surv_days,
          min_surv_people = min_surv_people,
          output_dir = output_dir,
          output_name = paste0(mutation_type,"_piu_cdr_univariate.tsv"))
        
      }else{
        survival_model_no_gender(
          surv_info_data = piu_info,
          interest_variable_info = piu_unite[[3]],
          min_surv_time = min_surv_days,
          min_surv_people = min_surv_people,
          output_dir = output_dir,
          output_name = paste0(mutation_type,"_piu_cdr_univariate.tsv"))
      }
    }else{
      cat("...this type of PIU level data not available.", "\n")
      
    }
    
    
  }
  
}



univariate_cox_model_for_lu = function(lu_filename,
                                                 cdr_clinical,
                                                 gender_as_covariate = T,
                                                 race_group_min = 6,
                                                 min_surv_days = 90,
                                                 min_surv_people = 5,
                                                 patient_sum_min = 3,
                                                 mutation_type = "somatic",
                                                 output_dir)
{
  
  
  
  if(file.exists(lu_filename))
  {
    

    
    lu_unite = gene_counts_cdr_clinical_unite(
      gene_count_filename = lu_filename,
      cdr_clinical = cdr_clinical,
      patient_sum_min = patient_sum_min,
      output_dir = output_dir,
      output_name = paste0(mutation_type,"_lu_cdr_clinical_unite.tsv"))
    
    
    
    if(nrow(lu_unite[[1]])>0)
    {
      lu_info = cdr_tidy_up_for_model(
        interest_variable_info = lu_unite[[2]],
        unite_data = lu_unite[[1]],
        race_group_min = race_group_min,
        output_dir = output_dir,
        output_name = paste0(mutation_type,"_lu_survival_info.tsv"))
      

      
      if(gender_as_covariate == T)
      {
        survival_model(
          surv_info_data = lu_info,
          interest_variable_info = lu_unite[[2]],
          min_surv_time = min_surv_days,
          min_surv_people = min_surv_people,
          output_dir = output_dir,
          output_name = paste0(mutation_type,"_lu_cdr_univariate.tsv"))
        
      }else{
        survival_model_no_gender(
          surv_info_data = lu_info,
          interest_variable_info = lu_unite[[2]],
          min_surv_time = min_surv_days,
          min_surv_people = min_surv_people,
          output_dir = output_dir,
          output_name = paste0(mutation_type,"_lu_cdr_univariate.tsv"))
      }
    }else{
      cat("...this lu data not available.", "\n")
      
    }
    

  }else{

  }
  

}





univariate_cox_model_for_ncu = function(ncu_filename,
                                                cdr_clinical,
                                                gender_as_covariate = T,
                                                race_group_min = 6,
                                                min_surv_days = 90,
                                                min_surv_people = 5,
                                                patient_sum_min = 3,
                                                mutation_type = "somatic",
                                                output_dir)
{
  
  if(file.exists(ncu_filename))
  {
    
    
    ncu_unite = gene_counts_cdr_clinical_unite(
      gene_count_filename = ncu_filename,
      cdr_clinical = cdr_clinical,
      patient_sum_min = patient_sum_min,
      output_dir = output_dir,
      output_name = paste0(mutation_type,"_ncu_cdr_clinical_unite.tsv"))
    
    
    
    if(nrow(ncu_unite[[1]])>0)
    {
      ncu_info = cdr_tidy_up_for_model(
        interest_variable_info = ncu_unite[[2]],
        unite_data = ncu_unite[[1]],
        race_group_min = race_group_min,
        output_dir = output_dir,
        output_name = paste0(mutation_type,"_ncu_survival_info.tsv"))
      

      if(length(ncu_unite[[2]])>0)
      {
        
        if(gender_as_covariate == T)
        {
          survival_model(
            surv_info_data = ncu_info,
            interest_variable_info = ncu_unite[[2]],
            min_surv_time = min_surv_days,
            min_surv_people = min_surv_people,
            output_dir = output_dir,
            output_name = paste0(mutation_type,"_ncu_cdr_univariate.tsv"))
          
        }else{
          survival_model_no_gender(
            surv_info_data = ncu_info,
            interest_variable_info = ncu_unite[[2]],
            min_surv_time = min_surv_days,
            min_surv_people = min_surv_people,
            output_dir = output_dir,
            output_name = paste0(mutation_type,"_ncu_cdr_univariate.tsv"))
        }
        
      }
      
      
    }else{
      cat("...this ncu data not available.", "\n")
      
    }
    


  }else{

  }
  

}


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
  
  # 
  # dss_data = unite_data %>%
  #   dplyr::select(barcode, race, DSS, DSS.time) %>%
  #   na.omit()
  # dss_race_freq = as.data.frame(table(dss_data$race))
  # dss_race_minor = dss_race_freq %>%
  #   dplyr::filter(Freq<race_group_min)
  # dss_data$race[which(dss_data$race %in% dss_race_minor$Var1)] = "OTHER"
  # colnames(dss_data)[which(colnames(dss_data)=="race")] = "dss_race"
  # 
  # 
  # dfi_data = unite_data %>%
  #   dplyr::select(barcode, race, DFI, DFI.time) %>%
  #   na.omit()
  # 
  # dfi_race_freq = as.data.frame(table(dfi_data$race))
  # dfi_race_minor = dfi_race_freq %>%
  #   dplyr::filter(Freq<race_group_min)
  # dfi_data$race[which(dfi_data$race %in% dfi_race_minor$Var1)] = "OTHER"
  # colnames(dfi_data)[which(colnames(dfi_data)=="race")] = "dfi_race"
  # 
  # 
  # pfi_data = unite_data %>%
  #   dplyr::select(barcode, race, PFI, PFI.time) %>%
  #   na.omit()
  # pfi_race_freq = as.data.frame(table(pfi_data$race))
  # pfi_race_minor = pfi_race_freq %>%
  #   dplyr::filter(Freq<race_group_min)
  # pfi_data$race[which(pfi_data$race %in% pfi_race_minor$Var1)] = "OTHER"
  # colnames(pfi_data)[which(colnames(pfi_data)=="race")] = "pfi_race"
  # 
  # 
  
  
  # patient_count_survival = patient_count %>%
  #   dplyr::left_join(os_data, by = "barcode") %>%
  #   dplyr::left_join(dss_data, by = "barcode") %>%
  #   dplyr::left_join(dfi_data, by =  "barcode") %>%
  #   dplyr::left_join(pfi_data, by = "barcode") %>%
  #   dplyr::select(barcode, age, gender,
  #                 os_race, OS, OS.time,
  #                 dss_race, DSS, DSS.time,
  #                 dfi_race, DFI, DFI.time,
  #                 pfi_race, PFI, PFI.time,
  #                 one_of(interest_variable_info))
  # 

  patient_count_survival = patient_count %>%
    dplyr::left_join(os_data, by = "Tumor_Sample_Barcode") %>%
      dplyr::select(Tumor_Sample_Barcode, age, gender,
                  os_race, OS, OS.time,
                  one_of(interest_variable_info))
  
  
  
  
  
 # write.table(patient_count_survival, paste0(output_dir, output_name),
   #           quote = F, row.names = F, sep = "\t")
  
  return(patient_count_survival)
  
}



piu_counts_cdr_clinical_unite = function(piu_count_filename,
                                         cdr_clinical,
                                         patient_sum_min,
                                         output_dir,
                                         output_name)
  
{
  
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
    dplyr::select(uniprot_accession, start_position, end_position, unit_name, gene_name, gene_id, unit_label,gene_info, piu_info, row_sum,
                  everything())
  

  if(nrow(piu_count_sel)>0)
  {
    
    col_seq = seq(1:ncol(piu_count_sel))
    which_barcode = col_seq[-c(1:10)]
    
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






survival_model= function(surv_info_data,
                                     interest_variable_info,
                                     min_surv_time,
                                     min_surv_people,
                                     output_dir,
                                     output_name)
{
  
  endpoint_flag = c(T)
  v_status = c("OS")
  v_time = c("OS.time")
  v_race = c("os_race")
  
  
  q_status = quos(OS)
  q_time = quos(OS.time)
  q_race = quos(os_race)
  
  
  list_surv_models = vector(mode = "list", length = 1)
  
  for(i in 1:1)
  {
    
   # i = 1
    this_race = q_race[[i]]
    this_status = q_status[[i]]
    this_time = q_time[[i]]
    
    
    this_race_v = v_race[i]
    this_status_v = v_status[i]
    this_time_v = v_time[i]
    
    
    
    this_surv_data = surv_info_data %>%
      dplyr::select(Tumor_Sample_Barcode, age,gender, !!this_race, !!this_status, !!this_time) %>%
      na.omit()%>%
      dplyr::arrange(!!this_status, !!this_time) %>%
      dplyr::filter(!!this_time >= min_surv_time)
    
    if(nrow(this_surv_data) == 0)
    {
      endpoint_flag[i] = F
    }else{
      
      this_surv_table = as.data.frame(table(this_surv_data[[this_status_v]]))
      
      if(min(this_surv_table$Freq)< min_surv_people)
      {
        endpoint_flag[i] = F
      }else{
        this_age_data = as.numeric(this_surv_data$age)
        this_gender_data = relevel(as.factor(this_surv_data$gender), ref = unique(this_surv_data$gender)[1])
        this_race_data = relevel(as.factor(this_surv_data[[this_race_v]]), ref = unique(this_surv_data[[this_race_v]])[1])
        this_surv_object = Surv(time = this_surv_data[[this_time_v]], event = this_surv_data[[this_status_v]])
        
      }
    }
    
    
    this_surv_result_df = rbindlist(lapply(1:length(interest_variable_info), function(x)
    {
      this_count_df = surv_info_data %>%
        dplyr::select(Tumor_Sample_Barcode, one_of(interest_variable_info[x]))
      this_count = this_surv_data %>%
        dplyr::left_join(this_count_df, by = "Tumor_Sample_Barcode")
      this_num_patients = sum(this_count[,7]!= 0)
      this_total_patients = nrow(this_count)
      this_count_coeff = NA
      this_count_exp_coeff = NA
      this_count_pval = NA
      this_assum_test = NA
      
      if(endpoint_flag[i] == T)
      {
        if(length(unique(this_gender_data))>1)
        {
          if(length(unique(this_race_data))>1)
          {
            this_model = coxph(this_surv_object ~  this_age_data  + this_gender_data + this_race_data + 
                                 this_count[,7])   
          }else{
            this_model = coxph(this_surv_object ~  this_age_data  + this_gender_data +
                                 this_count[,7])   
          }
          
        }else{
          if(length(unique(this_race_data))>1)
          {
            this_model = coxph(this_surv_object ~  this_age_data  + this_race_data + 
                                 this_count[,7])   
          }else{
            this_model = coxph(this_surv_object ~  this_age_data +
                                 this_count[,7])   
          }
          
        }
        
        
        this_test = cox.zph(this_model)
        this_table = this_test$table
        
        
        
        
        this = summary(this_model)
        this_coef = this$coefficients
        r_this =  length(unique(this_gender_data)) +length(unique(this_race_data))
        this_count_coeff = this_coef[r_this,1]
        this_count_exp_coeff = this_coef[r_this,2]
        this_count_pval = this_coef[r_this,5]
        this_assum_test = this_table[r_this,3]
        
        
        
      }
      
      one_surv_result = data.frame(count_info = interest_variable_info[x],
                                   this_num_patients, this_total_patients,
                                   this_count_coeff, this_count_exp_coeff, this_count_pval,
                                   this_assum_test,
                                   this_count_qval = NA,
                                   stringsAsFactors = F)
      
      if(x%%100 == 0)
        cat(x, i, "\n")
      
      return(one_surv_result)
      
    }))
    
    
    
    if(length(interest_variable_info)>5)
    {
      if(endpoint_flag[i] == T)
      {
        if(length(this_surv_result_df$this_count_pval)>300 & min(this_surv_result_df$this_count_pval,na.rm = T)<0.05 & max(this_surv_result_df$this_count_pval,na.rm = T)>0.95)
        {
          this_qval = qvalue(this_surv_result_df$this_count_pval)
        }else{
          this_qval = qvalue(this_surv_result_df$this_count_pval, pi0 = 1)
          
        }
        this_surv_result_df$this_count_qval = this_qval$qvalues
        
      }
    }else{
      cat("Small unit size, no q-val estimation.", "\n")
      
    }
    
    
    old_colnames = colnames(this_surv_result_df)
    new_colnames = gsub("this", v_status[i], old_colnames)
    colnames(this_surv_result_df) = new_colnames
    
    
    list_surv_models[[i]] = this_surv_result_df
    
    
  }
  
  
  
  final_result = list_surv_models %>%
    Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="count_info"), .) %>%
    dplyr::arrange(desc(OS_num_patients)) %>%
    replace(is.na(.),"")
  
  
  write.table(final_result, paste0(output_dir, output_name),
              quote = F, row.names = F, sep = "\t")
  
  
}









survival_model_no_gender= function(surv_info_data,
                                               interest_variable_info,
                                               min_surv_time,
                                               min_surv_people,
                                               output_dir,
                                               output_name)
{

  
  endpoint_flag = c(T)
  v_status = c("OS")
  v_time = c("OS.time")
  v_race = c("os_race")
  
  q_status = quos(OS)
  q_time = quos(OS.time)
  q_race = quos(os_race)
  
  
  list_surv_models = vector(mode = "list", length = 1)
  
  for(i in 1:1)
  {

    this_race = q_race[[i]]
    this_status = q_status[[i]]
    this_time = q_time[[i]]
    
    
    this_race_v = v_race[i]
    this_status_v = v_status[i]
    this_time_v = v_time[i]
    
    
    
    this_surv_data = surv_info_data %>%
      dplyr::select(Tumor_Sample_Barcode, age,gender, !!this_race, !!this_status, !!this_time) %>%
      na.omit()%>%
      dplyr::arrange(!!this_status, !!this_time) %>%
      dplyr::filter(!!this_time >= min_surv_time)
    
    if(nrow(this_surv_data) == 0)
    {
      endpoint_flag[i] = F
    }else{
      
      this_surv_table = as.data.frame(table(this_surv_data[[this_status_v]]))
      
      if(min(this_surv_table$Freq)< min_surv_people)
      {
        endpoint_flag[i] = F
      }else{
        this_age_data = as.numeric(this_surv_data$age)
        this_gender_data = relevel(as.factor(this_surv_data$gender), ref = unique(this_surv_data$gender)[1])
        this_race_data = relevel(as.factor(this_surv_data[[this_race_v]]), ref = unique(this_surv_data[[this_race_v]])[1])
        this_surv_object = Surv(time = this_surv_data[[this_time_v]], event = this_surv_data[[this_status_v]])
        
      }
    }
    
    
    this_surv_result_df = rbindlist(lapply(1:length(interest_variable_info), function(x)
    {

      this_count_df = surv_info_data %>%
        dplyr::select(Tumor_Sample_Barcode, one_of(interest_variable_info[x]))
      this_count = this_surv_data %>%
        dplyr::left_join(this_count_df, by = "Tumor_Sample_Barcode")
      this_num_patients = sum(this_count[,7]!= 0)
      this_total_patients = nrow(this_count)
      this_count_coeff = NA
      this_count_exp_coeff = NA
      this_count_pval = NA
      this_assum_test = NA
      
      
      if(endpoint_flag[i] == T)
      {
        
        if(length(unique(this_race_data))>1)
        {
          this_model = coxph(this_surv_object ~  this_age_data + this_race_data + 
                               this_count[,7])   
          
          
        }else{
          this_model = coxph(this_surv_object ~  this_age_data +
                               this_count[,7])   
        }
        
        
        this_test = cox.zph(this_model)
        this_table = this_test$table
        
        this = summary(this_model)
        this_coef = this$coefficients
        r_this =  1 +length(unique(this_race_data))
        this_count_coeff = this_coef[r_this,1]
        this_count_exp_coeff = this_coef[r_this,2]
        this_count_pval = this_coef[r_this,5]
        this_assum_test = this_table[r_this,3]
        
        
      }
      
      one_surv_result = data.frame(count_info = interest_variable_info[x],
                                   this_num_patients, this_total_patients,
                                   this_count_coeff, this_count_exp_coeff, this_count_pval,
                                   this_assum_test,
                                   this_count_qval = NA,
                                   stringsAsFactors = F)
      
      if(x%%100 == 0)
        cat(x, i, "\n")
      
      return(one_surv_result)
      
    }))
    
    
    
    if(length(interest_variable_info)>5)
    {
      if(endpoint_flag[i] == T)
      {
        if(length(this_surv_result_df$this_count_pval)>300 & min(this_surv_result_df$this_count_pval,na.rm = T)<0.05 & max(this_surv_result_df$this_count_pval,na.rm = T)>0.95)
        {
          this_qval = qvalue(this_surv_result_df$this_count_pval)
        }else{
          this_qval = qvalue(this_surv_result_df$this_count_pval, pi0 = 1)
          
        }
        this_surv_result_df$this_count_qval = this_qval$qvalues
        
      }
    }else{
      
      cat("Small unit size, no q-val estimation.", "\n")
      
    }
    
    
    old_colnames = colnames(this_surv_result_df)
    new_colnames = gsub("this", v_status[i], old_colnames)
    colnames(this_surv_result_df) = new_colnames
    
    
    list_surv_models[[i]] = this_surv_result_df
    
    
  }
  
  
  
  final_result = list_surv_models %>%
    Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="count_info"), .) %>%
    dplyr::arrange(desc(OS_num_patients)) %>%
    replace(is.na(.),"")
  
  
  write.table(final_result, paste0(output_dir, output_name),
              quote = F, row.names = F, sep = "\t")
  
  
}






univariate_survival_significance_adjust_totalMut_stage = function(filter_surv_data,
                                                                  output_dir,
                                                                  surv_name,
                                                                  data_fold)
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
  
  
  
  output_name = paste0(surv_name, "_", data_fold, ".tsv")
  
  return(surv_result_df) 
  
}



