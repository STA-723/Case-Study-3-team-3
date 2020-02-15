DATA LIST FILE= "C:/Users/Frances/Documents/First_Year_Classes/Second/Case Studies/Case-Study-3-team-3/mice_imp1/mice_imp1__SPSS.txt"  free (TAB)
   / Imputation_ C1_A C1_B C2 C3_A C3_B C3_C C3_D C4 C5 C6 
  C7 C8 C9 C10_A C10_B C10_C C10_D C10_E C10_F C10_G C10_H 
  C10_I C10_J C10_K C11_A C11_B C11_C C11_D C11_E C12_A 
  C12_B C12_C C13_A C13_B C13_C C13_D C14 C15 C16_A C16_B 
  C16_C C16_D C16_E C16_F C16_G C16_H C16_I C16_J C16_K 
  C17_A C17_B C17_C C17_D C17_E C17_F C17_G C17_H C17_I 
  C17_J C17_K C17_L C18_A C18_B C18_C C18_D C18_E C18_F 
  C18_G C18_H C18_I C18_J C19_A C19_B C19_C C19_D C19_E 
  C19_F C20 C21 C22_A C22_B C22_C C22_D C22_E C22_F C22_G 
  C22_H C22_I C22_J C22_K C22_L C22_M C22_N .


VARIABLE LABELS
  Imputation_ "Imputation_" 
 C1_A "C1_A" 
 C1_B "C1_B" 
 C2 "C2" 
 C3_A "C3_A" 
 C3_B "C3_B" 
 C3_C "C3_C" 
 C3_D "C3_D" 
 C4 "C4" 
 C5 "C5" 
 C6 "C6" 
 C7 "C7" 
 C8 "C8" 
 C9 "C9" 
 C10_A "C10_A" 
 C10_B "C10_B" 
 C10_C "C10_C" 
 C10_D "C10_D" 
 C10_E "C10_E" 
 C10_F "C10_F" 
 C10_G "C10_G" 
 C10_H "C10_H" 
 C10_I "C10_I" 
 C10_J "C10_J" 
 C10_K "C10_K" 
 C11_A "C11_A" 
 C11_B "C11_B" 
 C11_C "C11_C" 
 C11_D "C11_D" 
 C11_E "C11_E" 
 C12_A "C12_A" 
 C12_B "C12_B" 
 C12_C "C12_C" 
 C13_A "C13_A" 
 C13_B "C13_B" 
 C13_C "C13_C" 
 C13_D "C13_D" 
 C14 "C14" 
 C15 "C15" 
 C16_A "C16_A" 
 C16_B "C16_B" 
 C16_C "C16_C" 
 C16_D "C16_D" 
 C16_E "C16_E" 
 C16_F "C16_F" 
 C16_G "C16_G" 
 C16_H "C16_H" 
 C16_I "C16_I" 
 C16_J "C16_J" 
 C16_K "C16_K" 
 C17_A "C17_A" 
 C17_B "C17_B" 
 C17_C "C17_C" 
 C17_D "C17_D" 
 C17_E "C17_E" 
 C17_F "C17_F" 
 C17_G "C17_G" 
 C17_H "C17_H" 
 C17_I "C17_I" 
 C17_J "C17_J" 
 C17_K "C17_K" 
 C17_L "C17_L" 
 C18_A "C18_A" 
 C18_B "C18_B" 
 C18_C "C18_C" 
 C18_D "C18_D" 
 C18_E "C18_E" 
 C18_F "C18_F" 
 C18_G "C18_G" 
 C18_H "C18_H" 
 C18_I "C18_I" 
 C18_J "C18_J" 
 C19_A "C19_A" 
 C19_B "C19_B" 
 C19_C "C19_C" 
 C19_D "C19_D" 
 C19_E "C19_E" 
 C19_F "C19_F" 
 C20 "C20" 
 C21 "C21" 
 C22_A "C22_A" 
 C22_B "C22_B" 
 C22_C "C22_C" 
 C22_D "C22_D" 
 C22_E "C22_E" 
 C22_F "C22_F" 
 C22_G "C22_G" 
 C22_H "C22_H" 
 C22_I "C22_I" 
 C22_J "C22_J" 
 C22_K "C22_K" 
 C22_L "C22_L" 
 C22_M "C22_M" 
 C22_N "C22_N" 
 .

EXECUTE.
SORT CASES by Imputation_.
SPLIT FILE layered by Imputation_.