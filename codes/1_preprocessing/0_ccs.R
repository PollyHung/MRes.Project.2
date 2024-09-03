## Index 
# 1. SAMPLE_LIST_PREPARE: Prepare sample lists for ccs and subsequent runs. 






## prepare sample list: ========================================================
if(SAMPLE_LIST_PREPARE){
  sample_list <- list.files("/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw/")
  text_files <- !grepl(".txt", sample_list) ## which files are NOT text files? 
  sample_list <- sample_list[text_files]
  write.table(sample_list, "codes/1_preprocessing/sample_list.txt", 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  sample_list <- setdiff(sample_list, "0700055A") ## exlude this sample because we have already run ccs on it
  sets <- split(sample_list, ceiling(seq_along(sample_list) / 4))
  for (i in 1:6) {
    cat(paste("Set", i, ":", paste(sets[[i]], collapse = ", "), "\n"))
    write.table(x = sets[[i]], file = paste0("codes/1_preprocessing/", paste0("Set_", i, ".txt")), 
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}






