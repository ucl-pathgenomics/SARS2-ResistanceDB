 takes the current db and returns a yaml file for the variant pipeline

# read the db
df = read.delim("rna_antiviral_db.tsv",sep = "\t")


# filter by mutation category
regex = "1-proven sars antiviral resistance|2-homologous nuc' resmut observed more than once"
keep = grepl(regex , df$judy_class)
df = df[keep,]


# filter by drug target
regex = "NSP12"
keep = grepl(regex , df$sarscov2_protein)
df = df[keep,]


# filter co-mutations
regex = "\\+"
keep = !grepl(regex , df$sarscov2_protein_pos)
df = df[keep,]


# get the unique mutations relative to orf1ab
unique_mutations = toupper(unique(df$sarscov2_protein_pos))


# generate yaml
text = readLines("ref/barebone.yaml")
insert = "      - "
AAs = c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")

for(mut in unique_mutations){
  
  
  if( grepl("[A-Z]{1}[0-9]{1,5}[A-Z]", mut) ){
    
    # if is not inferred
    out = paste0(insert , mut)
    text = c(text, out)
    next
    
  }else{
    
    # if is inferred, allow any Mutant Type, so insert them all
    # extract wt
    wt = stringr::str_extract(mut, "[A-Z]")
    
    # for all other amino acids
    for(AA in AAs[!grepl(wt , AAs)]){
      out = paste0(insert , paste0(mut, AA) )
      text = c(text, out)
    }
    
  }
  
}

# write the yaml file
writeLines(text, paste0("ref/dbForPipeline-", Sys.Date(), ".yaml"))
