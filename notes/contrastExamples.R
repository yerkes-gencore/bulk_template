# Example contrasts
# Dataset: longitudinal study with four groups and sex as a variable
# combined "group" and "timepoint" to one factor, looked for interactions with sex
bulk$md$design <- stats::model.matrix(~ 0 + grp.tPt * sex, data = bulk$dge$samples, contrasts.arg = list(sex = "contr.sum"))
colnames(bulk$md$design) <- gsub("grp.tPt", "", colnames(bulk$md$design)) # removing "grp.tPt" from names for brevity 
          ## (I am not doing the same with "sex" because sex is binary, and only one column "sex1" exists. 
                ## Make sure you know which sex is coded as -1 and 1)
colnames(bulk$md$design) <- gsub(":", ".i.", colnames(bulk$md$design)) # if you use interaction terms, replace colon with ".i." to avoid errors

contrastFunction <- function(x) {
  makeContrasts(
  # Group 1 difference from Day 0 
    g1_D001vsD000 = g1.D001 - g1.D000,
    g1_D056vsD000 = g1.D056 - g1.D000,
    g1_D057vsD000 = g1.D057 - g1.D000,
    
  # Group 1 difference from Day 56 
    g1_D057vsD056 = g1.D057 - g1.D056,
    
  # Group 1 sex differences
    # g1_D000_MvF = g1.D000.i.sex1, # Note: this is the intercept, so not a possible contrast
    g1_D001_MvF = g1.D001.i.sex1,
    g1_D056_MvF = g1.D056.i.sex1,

  # Difference between groups at specific timepoints
    D001_g1vsg2 = g1.D001 - g2.D001,
    D056_g1vsg2 = g1.D056 - g2.D056,
    
  # Difference between Group 1 and two control groups (g3 and g4)
    D001_g1vsg3g4 = g1.D001 - (g3.D001 + g4.D001)/2,
    
  # Difference of differences (how did the longitudinal changes differ between groups)
    D001vD000_g1vsg2 = (g1.D001 - g1.D000) - (g2.D001 - g2.D000),
    D0057vD056_g1vsg2 = (g1.D057 - g1.D056) - (g2.D057 - g2.D056),
    
  # Difference of sex-based interactions between groups
    MvF_D057_g1vsg2 = g1.D057.i.sex1 - g2.D057.i.sex1,

    levels = bulk$md$design) # right hand values must be included as colnames in the design!
}

## eval() then actually starts the function, running makeContrasts() as normal
bulk$md$contr.matrix <- eval(contrastFunction())

## Using the custom function getSourceText, we can print out the contrasts!
cat(getSourceText(contrastFunction,2,2)) # removes two lines from beginning and end
