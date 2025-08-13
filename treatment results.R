library(dplyr)
library(tidyr)

## 1) combined_all 
combined_all <- as.data.frame(t(as.matrix(cbind(beta_1_sub, beta_2_sub))))
# 
# rownames(combined_all) 

## 2)
sentrix_col <- if ("Sentrix_ID" %in% names(mars)) "Sentrix_ID" else "SentrixID"
mars[[sentrix_col]] <- as.character(mars[[sentrix_col]])

ids <- rownames(combined_all)
ids <- ids[!is.na(ids) & ids != ""]

mars_sub <- mars[mars[[sentrix_col]] %in% ids, , drop = FALSE]

cols  <- paste0("HMDme_0", 0:8)
weeks <- 0:8

subset <- mars_sub %>%
  select(all_of(sentrix_col), all_of(cols))

mars_sub_1 <- subset %>%
  mutate(across(all_of(cols), ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(across(all_of(cols), ~ replace_na(., 0)))

mars_sub_1$response_slope <- apply(
  mars_sub_1[, cols, drop = FALSE],
  1, function(x) coef(lm(x ~ weeks))[2]
)

# 
response_df <- mars_sub_1 %>%
  select(all_of(sentrix_col), response_slope) %>% distinct()
rownames(response_df) <- response_df[[sentrix_col]]

## 3)
analysis_all    <- list()
methylation_all <- list()

drug_names <- c("TCA","SSRI","SNRI","NASSA","NARI","SSRE","OTH","NL","PP","LT","BZD","SLP")

# 
rownames(mars_sub) <- as.character(mars_sub[[sentrix_col]])

for (r in drug_names) {
  
  # 
  drug_cols <- grep(r, colnames(mars_sub), value = TRUE)
  if (length(drug_cols) == 0) next    
  
  # 
  med <- mars_sub %>%
    select(all_of(drug_cols)) %>%
    mutate(
      NID  = mars_sub$NID,
      resp = mars_sub$SKRSME0W,  
      week = mars_sub$Sadweeks   
    )
  
  # 
  med[is.na(med)] <- 0
  
  # 
  med$med_in <- rowSums(med[, drug_cols, drop = FALSE], na.rm = TRUE)
  
  # 
  last3 <- tail(drug_cols, min(3, length(drug_cols)))
  
  filtered_med <- med[
    rowSums(med[, last3,     drop = FALSE], na.rm = TRUE) >= 1 &
      rowSums(med[, drug_cols, drop = FALSE], na.rm = TRUE) >= 4,
  ]
  
  med_nam <- rownames(filtered_med)
  if (length(med_nam) == 0) next
  
  # 
  response_med <- response_df[response_df[[sentrix_col]] %in% med_nam, , drop = FALSE]
  if (nrow(response_med) == 0) next
  
  #
  methylation_all[[r]] <- combined_all[
    rownames(combined_all) %in% rownames(response_med),
    , drop = FALSE
  ]
  
  analysis_all[[r]] <- response_med[
    rownames(response_med) %in% rownames(methylation_all[[r]]),
    , drop = FALSE
  ]
}
