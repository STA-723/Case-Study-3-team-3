# Library
library(dplyr)
library(tidyr)
library(mice, lib.loc = "~/R/x86_64-pc-linux-gnu-library/3.6")
library(miceadds, lib.loc = "~/R/x86_64-pc-linux-gnu-library/3.6")

# Load in data
read_cas <- function(folder, file, recordfile) {
  skip <- grep("--------", readLines(unz(folder, recordfile)))
  record_layout <- read.table(unz(folder, recordfile), 
                              fill = TRUE, skip = skip, header = FALSE)
  colnames(record_layout) <- c("variable_name", "start_col", "end_col", "type")
  record_layout <- record_layout %>% filter(!is.na(end_col))
  widths <- as.numeric(as.character(record_layout$end_col)) - as.numeric(as.character(record_layout$start_col)) + 1
  cas <- read.fwf(unz(folder, file), widths = widths, header	= FALSE,
                  col.names = record_layout$variable_name)
  return(cas)
}
cas97 <- read_cas(folder = "Harvard_CAS_1997.zip", 
                  file = "Harvard_CAS_1997/DS0001/03163-0001-Data.txt",
                  recordfile = "Harvard_CAS_1997/DS0001/03163-0001-Record_layout.txt")
print("Read in file")

# Fill in NA due to skipped questions and coding error
# Section A
cas97 <- cas97 %>%
  mutate(A8_answered = cas97 %>%
           select(which(grepl("^A8_", colnames(.)))) %>%
           rowSums(na.rm = TRUE)) %>%
  mutate(A8_1 = ifelse(is.na(A8_1) & A8_answered > 0, 0, A8_1),
         A8_2 = ifelse(is.na(A8_2) & A8_answered > 0, 0, A8_2),
         A8_3 = ifelse(is.na(A8_3) & A8_answered > 0, 0, A8_3),
         A8_4 = ifelse(is.na(A8_4) & A8_answered > 0, 0, A8_4)) %>%
  # A6, A7 should be dropped
  select(-A6, -A7) %>%
  # Drop A4, A5 that are transfer questions, otherwise have to dummify
  select(-A4, -A5, -A8_answered, -A8)

# Section B
cas97 <- cas97 %>%
  mutate(B10_answered = cas97 %>%
           select(which(grepl("^B10_", colnames(.)))) %>%
           rowSums(na.rm = TRUE)) %>%
  mutate(B10_1 = ifelse(is.na(B10_1) & B10_answered > 0, 0, B10_1),
         B10_2 = ifelse(is.na(B10_2) & B10_answered > 0, 0, B10_2),
         B10_3 = ifelse(is.na(B10_3) & B10_answered > 0, 0, B10_3),
         B10_4 = ifelse(is.na(B10_4) & B10_answered > 0, 0, B10_4),
         B10_5 = ifelse(is.na(B10_5) & B10_answered > 0, 0, B10_5),
         B10_6 = ifelse(is.na(B10_6) & B10_answered > 0, 0, B10_6),
         B10_7 = ifelse(is.na(B10_7) & B10_answered > 0, 0, B10_7),
         B10_8 = ifelse(is.na(B10_8) & B10_answered > 0, 0, B10_8)) %>%
  select(-B10_answered)

# Section C

# Section D
cas97 <- cas97 %>%
  select(-D7_1, -D7_2, -D7_3, -D7_4)

# Section E
cas97 <- cas97 %>%
  # Fill NA resulted from dummifying vars E23
  mutate(E23_answered = cas97 %>%
           select(which(grepl("^E23_", colnames(.)))) %>%
           rowSums(na.rm = TRUE)) %>%
  mutate(E23_1 = ifelse(is.na(E23_1) & E23_answered > 0, 0, E23_1),
         E23_2 = ifelse(is.na(E23_2) & E23_answered > 0, 0, E23_2),
         E23_3 = ifelse(is.na(E23_3) & E23_answered > 0, 0, E23_3)) %>%
  # Fill E27 with 0 for people less than 21 years old
  mutate(E27_A = ifelse(is.na(E27_A) & AGELT21 == 1, 0, E27_A)) %>%
  select(-E27_B, -E27_C) %>%
  # Combine E24, 25, 26: Did you get seriously injured within 6 hours of drinking: fill 0 in NA E25
  replace_na(list(E25 = 0)) %>%
  # E10-12 skipped if never had sex. Can combine E9 and E11, drop E10, E12
  mutate(E9n11 = ifelse(is.na(E11) & E9 == 1, 0, E11 + 1)) %>%
  select(-E9, -E10, -E12, -E11) %>%
  # E13-15 skipped if male
  replace_na(list(E13 = 0, E14 = 0, E15 = 0)) %>%
  # Drop E24, E26, E23_answered
  select(-E24, -E26, -E23_answered, -E23)

# Section F
cas97 <- cas97 %>%
  select(-F70FRND, -F30FRND)

# Section G
cas97 <- cas97 %>% 
  mutate(G3_answered = cas97 %>%
           select(which(grepl("^G3_", colnames(.)))) %>%
           rowSums(na.rm = TRUE)) %>%
  mutate(G3_1 = ifelse(is.na(G3_1) & G3_answered > 0, 0, G3_1),
         G3_2 = ifelse(is.na(G3_2) & G3_answered > 0, 0, G3_2),
         G3_3 = ifelse(is.na(G3_3) & G3_answered > 0, 0, G3_3),
         G3_4 = ifelse(is.na(G3_4) & G3_answered > 0, 0, G3_4),
         G3_5 = ifelse(is.na(G3_5) & G3_answered > 0, 0, G3_5),
  ) %>%
  select(-G3_answered, -G3) %>%
  mutate(G4_1 = I(G4==1),
         G4_2 = I(G4==2),
         G4_3 = I(G4==3),
         G4_4 = I(G4==4),
         G4_5 = I(G4==5),
         G4_6 = I(G4==6)) %>%  #indicators for religion
  select(-G4) %>%
  mutate(G14_NA = as.numeric(G14==5),
         G14_none = as.numeric(G14==6),
         G15_NA = as.numeric(G15==5),
         G15_none = as.numeric(G15==6)) %>% #indicators for don't know/NA, parental drinking
  mutate(G16_none=as.numeric(G16==4)) #indicator for no family agreement about drinking
print("Fill in some NAs")

# Select just the variables we use in the model
indep_vars_patterns <- "^B8_B$|^B8_C$|^B12_A$|^B7_A$|^B7_B$|^B7_C$|^B7_D$|^B7_E$|^B7_F$|^B8_D$|^B12_B$|^B9_C$|^B9_G$|^B1$|^B6_A$|^B6_B$|^B9_A$|^B9_D$|^B9_E$|^B9_F$|^B13_A$|^B13_B$|^B9_B$|^B11$|^B6_E$|^B6_C$|^B13_E$|^B2$|^B3$|^B13_C$|^B13_D$|^B6_D$|^F5_A|^F5_B|^F5_C|^F5_D|^F5_E|^F5_F|^F5_G|^F5_H|^F6_A|^F6_B|^F6_C|^F6_D|^F6_E|^F6_E|^F6_F|^F6_G|^F6_H|^F6_I|^F1|^F4|^F2|^F3|^G2|^G3|^G4|^G8|^G9|^G10|^G15|^G14|^G16|^G17|^G13|^A8|^D2|^D8|^D9|^D3|^D1|^A9|^A10"
cas97_small <- cas97 %>%
  select(which(grepl(indep_vars_patterns, colnames(.))))

# Run MICE
print("Run MICE")
mice_cas97 <- mice(cas97, m = 5)
write.mice.imputation(mi.res = mice_cas97, name = "mice_cas97" )
print("Finish!")




