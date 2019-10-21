


# a function to avoid file dir misform ------------------------------------

heal_dir = function(input_dir)
{
  len_input = nchar(input_dir)
  split_input = unlist(strsplit(input_dir, split = ""))
  
  if(split_input[len_input] != "/")
  {
    output_dir = paste0(input_dir, "/")
    
  }else{
    output_dir = input_dir
  }
  
  return(output_dir)
  
}





select_cancer_mc3 = function(mc3_df, cancer_barcode, output_dir, output_name)
  
  
{
   

  sel_mc3 = mc3_df %>%
    dplyr::filter(Tumor_Sample_Barcode %in% cancer_barcode) %>%
    unique()%>%
    dplyr::group_by(Hugo_Symbol, Gene, Chromosome, Start_Position, End_Position,
                    Variant_Classification, Variant_Type, HGVSc, HGVSp) %>%
    dplyr::summarise(agg_sample_id = paste(Tumor_Sample_Barcode, collapse = "_"), mut_freq = n()) %>%
    dplyr::arrange(desc(mut_freq))
  
  write.table(sel_mc3, paste0(output_dir, output_name),
              quote = F, row.names = F, sep = "\t")
  
  
}



# divide into pc and npc --------------------------------------------------

divide_somatic_to_pc_npc = function(mc3_data_name,
                                    output_dir,
                                    pc_output_name,
                                    npc_output_name)
{
  
 # mc3_data_name = "/data/ginny/tcga_pancan/important_files/mc3_info.tsv"
  
  mc3_data = fread(mc3_data_name, stringsAsFactors = F)
  
  
  ### here the goal is to find out the ones with interpretable HGVSp 
  
  
  pc_data = mc3_data %>%
    dplyr::filter(!grepl("Silent", Variant_Classification)) %>%
    dplyr::filter(nchar(HGVSp)>1) %>%
    dplyr::filter(! HGVSp == "p.=")## changed according to Dan's advise
  
  
  npc_data = mc3_data %>%
    dplyr::filter(!grepl("Silent", Variant_Classification)) %>%
    dplyr::filter(nchar(HGVSp)<=1)  ## changed according to Dan's advise
  
  ### no protein info is marked by a dot.
  
  write.table(pc_data, paste0(output_dir, pc_output_name), 
              quote = F, row.names = F, sep = "\t")
  
  
  write.table(npc_data, paste0(output_dir, npc_output_name), 
              quote = F, row.names = F, sep = "\t")
  
  
  
}





# annotate_position_for_pc ------------------------------------------------




annotate_mc3_pc_position_info= function(pc_data_name, 
                                        output_dir,
                                        output_name)
{
  pc_data = fread(pc_data_name, stringsAsFactors = F)
  
  
  p_p_pos = gsub("fs.*\\D*",  "",pc_data$HGVSp)
  
  
  p_pos =  gsub("[[:alpha:]]","",p_p_pos)
  fp_pos = gsub("\\.","",p_pos)
  
  start_end_pos = rbindlist(lapply(1:length(fp_pos), function(t){
    
    get_it = unlist(strsplit(fp_pos[t], split = "_"))
    
    start_end_df = data.frame(pos_start = 0, pos_end = 0)
    
    get_it[1] = gsub("\\D", "", get_it[1])
    
    
    start_end_df$pos_start = as.integer(get_it[1])
    
    
    if(length(get_it)!=2)
    {
      start_end_df$pos_end =  as.integer(get_it[1])
    }else{
      get_it[2] = gsub("\\D", "", get_it[2])
      start_end_df$pos_end =  as.integer(get_it[2])
    }
    
    return(start_end_df)
  }
  ))
  
  

  get_mc3 = pc_data%>%
    dplyr::mutate(prot_start_pos = start_end_pos$pos_start, 
                  prot_end_pos = start_end_pos$pos_end)%>%
    dplyr::select(Chromosome, Start_Position, End_Position, Variant_Classification, Variant_Type, 
                  Hugo_Symbol, Gene, HGVSc, HGVSp, prot_start_pos, prot_end_pos, agg_sample_id, mut_freq) %>%
    unique()
  
  
  
  write.table(get_mc3, paste0(output_dir, output_name),
              quote = F, row.names = F, sep = "\t")
  
  
  
}




# construct gene level count matrix --------------------------------------


gene_level_matrix = function(pc_data_name,
                              cancer_barcode,
                              mut_freq_min,
                              output_dir,
                              output_filename)
{
  # pc_data_name = pc_data_name
  # cancer_barcode = cancer_barcode
  # mut_freq_min = gene_level_mut_min
  # output_dir = output_dir
  # output_filename = "gene_level_count_matrix.tsv"
  # 
  
  pc_data = fread(pc_data_name, stringsAsFactors = F)
  
  
  gene_df = pc_data%>%
    dplyr::mutate(gene_info = paste(Hugo_Symbol,Gene, sep = "_"))
  
  
  all_gene = unique(gene_df$gene_info)
  
  sl_matrix = matrix(0, nrow = length(all_gene), ncol = length(cancer_barcode))
  colnames(sl_matrix) = cancer_barcode
  
  
  
  if(nrow(sl_matrix)>0)
  {
    for(i in 1:nrow(sl_matrix))
    {
      get_this_gene_df = gene_df%>%
        dplyr::filter(gene_info == all_gene[i])
      
      
      this_gene_patients = unlist(lapply(1:nrow(get_this_gene_df), function(x){
        this_patient =  unlist(strsplit(get_this_gene_df$agg_sample_id[x], split = "_"))
        return(this_patient)
      }))
      
      tdf = as.data.frame(table(this_gene_patients), stringsAsFactors = F)
      matrix_df = data.frame(pbr = cancer_barcode, count = 0, stringsAsFactors = F)
      
      get_counts = matrix_df %>%
        dplyr::left_join(tdf, by = c("pbr" = "this_gene_patients")) %>%
        replace(is.na(.),0) 
      

      sl_matrix[i,] = get_counts$Freq
      
      if(i%%100 ==0)
        cat(i,"\n")
      
    }
    
    
    
    gene_count_df = data.frame(all_gene, sl_matrix, row_sum = rowSums(sl_matrix), stringsAsFactors = F)
    colnames(gene_count_df) = c("gene_info", cancer_barcode, "row_sum")
    
    gene_count =  gene_count_df %>%
      dplyr::arrange(desc(row_sum)) %>%
      dplyr::filter(row_sum >= mut_freq_min)
    
    
    write.table(gene_count, paste0(output_dir, output_filename),
                quote =F, row.names = F, sep = "\t")
    
    
  }else{
    
    cat("No frequent gene level data.","\n")
  }
  
}







# construct locus level count matrix --------------------------------------


locus_level_matrix = function(pc_data_name,
                              cancer_barcode,
                              mut_freq_min,
                              output_dir,
                              output_filename)
{
   
    pc_data = fread(pc_data_name, stringsAsFactors = F)
  

  loci_df = pc_data%>%
    dplyr::mutate(loci_info = paste(Chromosome,Start_Position, End_Position, Hugo_Symbol,Gene,HGVSc, 
                                    HGVSp, prot_start_pos, prot_end_pos,sep = "_"))%>%
    dplyr::arrange(desc(mut_freq))  %>%
    dplyr::filter(mut_freq>=mut_freq_min)
  
  
  
  
  sl_matrix = matrix(0, nrow = nrow(loci_df), ncol = length(cancer_barcode))
  colnames(sl_matrix) = cancer_barcode
  
  if(nrow(sl_matrix)>0)
  {
    for(i in 1:nrow(sl_matrix))
    {
      get_sample = unlist(strsplit(loci_df$agg_sample_id[i], split = "_"))
      
      which_col = which(colnames(sl_matrix) %in% get_sample)
      
      sl_matrix[i,which_col] = 1
      
      if(i%%100 ==0)
        cat(i,"\n")
      
    }
    
    
    
    loci_count_df = data.frame(loci_df$loci_info, sl_matrix, row_sum = rowSums(sl_matrix), stringsAsFactors = F)
    colnames(loci_count_df) = c("locus_info", cancer_barcode, "row_sum")
    
    loci_count =  loci_count_df %>%
      dplyr::arrange(desc(row_sum))
    
    
    write.table(loci_count, paste0(output_dir, output_filename),
                quote =F, row.names = F, sep = "\t")
    
    
  }else{
    
    cat("No frequent locus level data.","\n")
  }
  
}




# map pc to piu and summarise to lu-------------------------------------



mc3_map_uni_piu = function(ptm_pfam_df,
                           pc_data_name,
                           cancer_barcode,
                           output_dir,
                           piu_output_filename,
                           lu_output_filename)
  
  
{
  

  pc_mut = fread(pc_data_name, stringsAsFactors = F)
  piu = ptm_pfam_df
  unique_prot = intersect(piu$gene_id, pc_mut$Gene)
  unique_sample = cancer_barcode
  
  watch_piu  = piu %>%
    dplyr::filter(gene_id %in% unique_prot)
  
  if(length(unique_prot) == 0)
  {
    cat("no single protien mapped", "\n")
  }else{
    
    
    piu_df = rbindlist(lapply(1:length(unique_prot), function(x)

    {
      
      this_prot = unique_prot[x]
      this_prot_piu = piu%>%
        dplyr::filter(gene_id == this_prot)

      count_piu_matrix = matrix(0, nrow = nrow(this_prot_piu), ncol = length(unique_sample))
      colnames(count_piu_matrix) = unique_sample
      
      local = pc_mut %>%
        dplyr::filter(Gene == this_prot)
      
      if(nrow(local)>0)
      {
        for(i in 1:nrow(local))
        {
          # i = 1
          get_sample = data.frame(this_seq = unlist(strsplit(local$agg_sample_id[i], split = "_")), stringsAsFactors = F)
          
          these_patient_col = which(unique_sample %in% get_sample$this_seq)
          
          piu_row_flag <- (local$prot_start_pos[i]>=this_prot_piu$start_position & local$prot_start_pos[i]<= this_prot_piu$end_position |local$prot_end_pos[i]>= this_prot_piu$start_position & local$prot_end_pos[i]<= this_prot_piu$end_position)
          
          
          count_piu_matrix[piu_row_flag, these_patient_col] =
            count_piu_matrix[piu_row_flag, these_patient_col] +1
          
          
        }
      }
      
      row_sum = rowSums(count_piu_matrix)
      
      
      this_prot_local_df = cbind(this_prot_piu, count_piu_matrix, row_sum)
      
      if(x%%100 == 0)
        cat(x, "\n")
      
      return(this_prot_local_df)
      
    }))
    
    with_info_piu_df = piu_df %>%
      dplyr::filter(row_sum>0) %>%
      dplyr::arrange(desc(row_sum))
    write.table(with_info_piu_df, paste0(output_dir, piu_output_filename), quote = F, row.names = F, sep = "\t")
    
    
    cat("PIU mapped","\n")
    

    lu_df = rbindlist(lapply(1:length(unique_prot), function(x)
      
    {
      
      this_prot = unique_prot[x]
      this_prot_piu = piu%>%
        dplyr::filter(gene_id == this_prot)
      
      this_lu_info = this_prot_piu %>%
        dplyr::select(gene_id,gene_name)%>%
        unique()%>%
        dplyr::top_n(n = 1, wt = gene_name)
      
      
      count_lu_matrix = matrix(0, nrow = 1, ncol = length(unique_sample))
      colnames(count_lu_matrix) = unique_sample
      
      local = pc_mut %>%
        dplyr::filter(Gene == this_prot)
      
      if(nrow(local)>0)
      {
        for(i in 1:nrow(local))
        {
          # i = 1
          get_sample = data.frame(this_seq = unlist(strsplit(local$agg_sample_id[i], split = "_")), stringsAsFactors = F)
          
          
          these_patient_col = which(unique_sample %in% get_sample$this_seq)
          
          
          piu_row_flag <- (local$prot_start_pos[i]>=this_prot_piu$start_position & local$prot_start_pos[i]<= this_prot_piu$end_position |local$prot_end_pos[i]>= this_prot_piu$start_position & local$prot_end_pos[i]<= this_prot_piu$end_position)
          
          
          if(sum(piu_row_flag)==0)
          {
            count_lu_matrix[1, these_patient_col] = count_lu_matrix[1, these_patient_col] +1
          }
          
          
        }
      }
      
      row_sum = rowSums(count_lu_matrix)
      
      
      this_prot_lu_df = cbind(this_lu_info, count_lu_matrix, row_sum)
      
      if(x%%100 == 0)
        cat(x, "\n")
      
      return(this_prot_lu_df)
      
    }))
    
    with_info_lu_df = lu_df %>%
      dplyr::filter(row_sum>0) %>%
      dplyr::arrange(desc(row_sum))
    write.table(with_info_lu_df, paste0(output_dir, lu_output_filename), quote = F, row.names = F, sep = "\t")
    
    
    
  }
  
}


# summarise for npc -------------------------------------------------------




mc3_map_npc = function(npc_data_name,
                       cancer_barcode,
                       output_dir,
                       output_filename)
  
  
{ 
  
  
  
  non_impact = fread(npc_data_name, stringsAsFactors = F)
  gene = non_impact %>%
    dplyr::select(Hugo_Symbol, Gene)%>%
    unique()
  
  unique_sample = cancer_barcode
  
  if(nrow(gene) == 0)
  {
    cat("no single protein mapped", "\n")
  }else{
    
    gene_df = rbindlist(lapply(1:nrow(gene), function(x)

    {
      this_gene_id = gene$Gene[x]
      
      
      count_gene = gene[x,]
      
      count_gene_matrix = matrix(0, nrow = 1, ncol = length(unique_sample))
      colnames(count_gene_matrix) = unique_sample
      
      this_non_impact = non_impact %>%
        dplyr::filter(Gene == this_gene_id)
      
      if(nrow(this_non_impact)>0)
      {
        for(i in 1:nrow(this_non_impact))
        {
          # i = 1
          get_sample = data.frame(this_seq = unlist(strsplit(this_non_impact$agg_sample_id[i], split = "_")), stringsAsFactors = F)
          
          these_patient_col = which(colnames(count_gene_matrix) %in% get_sample$this_seq)
          
          count_gene_matrix[1, these_patient_col] = 
            count_gene_matrix[1, these_patient_col] +1
          
        }
      }
      
      row_sum = rowSums(count_gene_matrix)
      
      
      this_gene_df = cbind(count_gene, count_gene_matrix, row_sum)
      
      if(x%%100 == 0)
        cat(x, "\n")
      
      return(this_gene_df)
      
    }))
    
    
    with_info_gene_df = gene_df %>%
      dplyr::filter(row_sum>0) %>%
      dplyr::arrange(desc(row_sum))%>%
      dplyr::rename(gene_name = Hugo_Symbol)%>%
      dplyr::rename(gene_id = Gene)
    
    
    write.table(with_info_gene_df, paste0(output_dir, output_filename), quote = F, row.names = F, sep = "\t")
    
    
  }
  
}




