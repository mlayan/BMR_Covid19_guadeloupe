periods = c(
  "Période pré-COVID (01/01/2020-31/01/2020)", 
  "3e vague COVID (15/04/2021-14/05/2021)", 
  "4e vague COVID (07/08/2021-06/09/2021)", 
  "Période inter-vague (01/05/2022-31/05/2022)"
  )

periods_start = as.Date(c(
  "2020-01-01", 
  "2021-04-15", 
  "2021-08-07", 
  "2022-05-01"
  ))

periods_end = as.Date(c(
  "2020-01-31", 
  "2021-05-14", 
  "2021-09-06", 
  "2022-05-31"
  ))

period_names_bold = c(
  "**Période pré-COVID\n(01/01/2020-31/01/2020)**", 
  "**3e vague COVID\n(15/04/2021-14/05/2021)**",
  "**4e vague COVID\n(07/08/2021-06/09/2021)**",
  "**Période inter-vague\n(01/05/2022-31/05/2022)**"
)

period_names = c(
  "Période pré-COVID\n(01/01/2020-31/01/2020)", 
  "3e vague COVID\n(15/04/2021-14/05/2021)",
  "4e vague COVID\n(07/08/2021-06/09/2021)",
  "Période inter-vague\n(01/05/2022-31/05/2022)"
)

period_names_short = c(
  "preCOVID", 
  "wave3",
  "wave4",
  "interwave"
)

bacteria_species = c(
  "PVTKP" = "K pneumoniae", 
  "PVTECO" = "E coli", 
  "PVTECC" = "E cloacae", 
  "PVTEAG" = "E aerogenes", 
  "PVTBA1" = "Other"
  )

sectors = c(
  "EXTENSION",
  "MAHOGANY",
  "HIBISCUS",
  "SELF",
  "FLAMBOYANT"
)

atb_class = c(
  "ACIDE CLAVULANIQUE" = "Penicillins",           
  "AMIKACINE" = "Aminosides",
  "AMIKLIN"  = "Aminosides",
  "AMINOSIDE" = "Aminosides",
  "AMOXICILIINE" = "Penicillins", 
  "AMOXICILLINE" = "Penicillins", 
  "AUGMENTIN" = "Penicillins", 
  "CEFAZOLINE" = "Cephalosporins",
  "CEFEPIME" = "Cephalosporins",
  "CEFOTAXIM" = "Cephalosporins",
  "CEFOTAXIME" = "Cephalosporins",
  "CEFTRIAXONE" = "Cephalosporins",
  "CEFUROXIME" = "Cephalosporins",
  "CLAFORAN" = "Cephalosporins",
  "CLINDAMYCINE" = "Macrolides",
  "FLAGYL" = "Nitroimidazole",
  "GENTAMICINE" = "Aminosides",
  "GENTAMYCINE" = "Aminosides",
  "LEVOFLOXACINE" = "Quinolones",
  "MEROPENEME" = "Cabapenems",
  "METRONIDAZOLE" = "Nitroimidazole",
  "OXACILLINE" = "Penicillins", 
  "PIPERACILLINE" = "Penicillins", 
  "RIFAMPICINE" = "Rifamycines",
  "ROCEPHINE" = "Cephalosporins",
  "ROVAMYCINE ROVAMYCINE ET LINEZOLID" = "Macrolides",
  "SPIRAMYCINE" = "Macrolides",
  "SULFAMETHOXAZOLE" = "Sulfamides",
  "TACIZOCILLINE" = "Penicillins",
  "TAZOBACTAM" = "Penicillins", 
  "TAZOCILLINE" = "Penicillins",
  "TRIMETHOPRIME" = "Sulfamides",
  "TRIMETOPRIME" = "Sulfamides",
  "VANCOMYCINE" = "Glycopeptides"
)

inf_loc_french = c(
  "0"= "Pvt rectal",
  "1"="Urine",
  "2"="Sang périphérique",
  "3"="Pus",
  "4"="Selles",
  "5"="LCR",
  "6"="Vaginal",
  "7"="Sperme",
  "8"="Gorge", 
  "9"="Crachat",
  "10"="LBA",
  "11"="Aspiration endo-trachéale",
  "12"="Liquide de ponction articulaire ou synoviale",
  "13"="Liquide pleural",
  "14"="Liquide péritonéal",
  "15"="Liquide d'ascite",
  "16"="Liquide péricardique",
  "17"="Cathéter",
  "18"="Drain",
  "19"="Sonde",
  "20"="Nasopharynx",
  "21"="PDP",
  "22"="Aspiration bronchique",
  "23"="Moelle osseuse",
  "24"="Aspiration duodénale",
  "25"="Prélèvement ORL",
  "26"="Sang cathéter"
)


inf_loc = c(
  "0"= "Rectal sample",
  "1"="Urine",
  "2"="Peripheral blood",
  "3"="Pus",
  "4"="Stools",
  "5"="Cerebrospinal fluid",
  "6"="Vaginal",
  "7"="Sperme",
  "8"="Throat", 
  "9"="Sputum",
  "10"="Bronchoalveolar lavage",
  "11"="Endo-tracheal aspiration",
  "12"="Joint or synovial puncture fluid",
  "13"="Pleural fluid",
  "14"="Peritoneal fluid",
  "15"="Ascites fluid",
  "16"="Pericardial fluid",
  "17"="Catheter",
  "18"="Drain",
  "19"="Line",
  "20"="Nasopharynx",
  "21"="Protected specimen brush",
  "22"="Bronchial aspiration",
  "23"="Bone marrow",
  "24"="Duodenal aspiration",
  "25"="ENT sample",
  "26"="Blood catheter" 
)

inf_loc_type = c(
  "0"= "Digestive",
  "1"="Urine",
  "2"="Blood",
  "3"="Skin",
  "4"="Digestive",
  "5"="Central nervous system",
  "6"="Genital",
  "7"="Genital",
  "8"="Respiratory", 
  "9"="Respiratory",
  "10"="Respiratory",
  "11"="Respiratory",
  "12"="Joint or synovial puncture fluid",
  "13"="Body fluids",
  "14"="Body fluids",
  "15"="Body fluids",
  "16"="Body fluids",
  "17"="Blood",
  "18"="Blood",
  "19"="Blood",
  "20"="Respiratory",
  "21"="Digestive",
  "22"="Respiratory",
  "23"="Central nervous system",
  "24"="Digestive",
  "25"="ENT",
  "26"="Blood" 
)