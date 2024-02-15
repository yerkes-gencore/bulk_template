## Cannot be sourced / ran as background due to need to restart R
if (!require('renv')) {
  install.packages('renv')
}
library(renv)

## Link for non-locked CRAN
## "https://cran.r-project.org" 
CRAN_mirror <- "https://packagemanager.posit.co/cran/__linux__/jammy/2023-10-30"
github_repos <- c('https://github.com/yerkes-gencore/gencoreBulk')
bioconductor <- TRUE ## you can also set this to a specific version

# repos <- BiocManager::repositories()
# repos["CRAN"] <- CRAN_mirror
# options(repos=repos)

init_settings <- list(
  snapshot.type = 'implicit', ## explicit = Use the description file to initialize the project
  use.cache = TRUE,
  ignored.packages = c()      ## Specify any packages to exclude from renv
)

renv::init(bioconductor = bioconductor, # repos = c(CRAN_mirror, github_repos), repos = repos
           settings = init_settings, bare = FALSE,
           load = TRUE, restart = TRUE, force = FALSE)
renv::activate()
renv::snapshot()
