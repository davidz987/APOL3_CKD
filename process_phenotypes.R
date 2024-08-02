### Code framework for generate phecode case/controls as well as QC for quantitative traits
### Patient-sensitive and biobank-protected data removed


library(dplyr)
library(tidyr)
library(readr)

###################
### Process phecodes
###################

### Read in ICD codes in long format
icd_codes <- data.table::fread("icd_codes_long.txt")

### Read in ICD (both 9 and 10) mappings to phecodes - from phewascatalog.org
icd_phecode_map <- read_delim("ICD-CM_phecode_unrolled.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

### Map ICD codes to phecodes
phecode9 <- icd_codes %>% 
        filter(code_standard_name == 'ICD9CM') %>% 
        select(ID, code, visit_id) %>% 
        unique() %>% 
        merge(icd_phecode_map %>% filter(flag == 9), by.x='code', by.y='ICD', allow.cartesian=TRUE) %>% 
        select(ID, phecode, visit_id) %>% 
        unique()

phecode10 <- icd_codes %>% 
        filter(code_standard_name == 'ICD10CM') %>% 
        select(ID, code, visit_id) %>% 
        unique() %>% 
        merge(icd_phecode_map %>% filter(flag == 10), by.x='code', by.y='ICD', allow.cartesian=TRUE) %>% 
        select(ID, phecode, visit_id) %>% 
        unique()

phecode_long <- phecode9 %>% rbind(phecode10)

phecode_mat <- phecode_long %>% 
        group_by(ID, phecode) %>% 
        summarise(num_visits = length(unique(visit_id)), .groups='drop') %>% 
        pivot_wider(names_from = phecode,
                values_from = num_visits,
                values_fill = 0)

phecode_mat <- phecode_mat %>% 
        select(c(ID, sort(colnames(phecode_mat[ , !colnames(phecode_mat) %in% c('ID')])))) %>% 
        arrange(ID)

### Generate case/controls using a rule of 2 - patient requires at least 2 instances of a diagnosis to be designated a case
filter_ro2 <- function(df) {
        print("Applying rule of 2...")
        df_R <- df[,-1]                         # Remove ID column to apply rule of 2 - add back later
        df_R[df_R == 1] <- NA
        df_R[df_R >= 2] <- 1
        df_R$ID <- df$ID
        return(df_R)
}

phecodes_R <- filter_ro2(phecode_mat) %>%
        select(c("ID"), everything())

write.table(phecodes_R, 
        file="Phecode_ro2.txt", 
        row.names=FALSE, quote=FALSE, sep="\t")


###################
### Process quantitative traits
###################

RankNorm <- function(col) {
        k = 0.375
        n <- length(col)
        r <- rank(col, ties.method = "average")
        out <- qnorm((r - k) / (n - (2 * k) + 1))
        return(out)
}

### Read in labs in long format
labs <- read.delim("labs_long.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)

### Apply filtering
### - Filter for Inpatient/Outpatient labs
### - Remove NAs, NULLs, infinite, non-numeric
### - Remove duplicates
### - Remove outliers (>4 SDs)
labs_formatted <- labs %>%
        filter(grepl("Outpatient", patient_class) | grepl("Inpatient", patient_class)) %>%
        select(ID, order_date, value) %>%
        mutate(across(where(is.character), ~na_if(., "NULL"))) %>%
        na.omit() %>%
        distinct() %>%
        filter(value >= mean(value) - (4 * sd(value))) %>%
        filter(value <= mean(value) + (4 * sd(value))) 

### - Calculate min, median, max
### - Inverse normal quantile transformation
labs_formatted <- labs_formatted %>%
        select(ID, value) %>%
        group_by(ID) %>%
        summarize_all(list(original_median = ~ median(., na.rm = TRUE), original_max = ~ max(., na.rm = TRUE), original_min = ~ min(., na.rm = TRUE))) %>%
        mutate(median = RankNorm(original_median)) %>%
        mutate(max = RankNorm(original_max)) %>%
        mutate(min = RankNorm(original_min))

write.table(labs_formatted, 
        file="qc_labs.txt", 
        row.names=FALSE, quote=FALSE, sep="\t")






