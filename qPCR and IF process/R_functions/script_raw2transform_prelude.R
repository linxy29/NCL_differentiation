knitr::opts_chunk$set(echo = TRUE, results = "markup", error = TRUE, fill = TRUE, out.width = "100%")
options(knitr.kable.NA = "", width = 1000)

load_packs  <- c("knitr","data.table", "DT", "stringr", "validate")
require_packs <- c("tools", "readxl")
import_packs <- c(load_packs, require_packs)

notinstall <- import_packs[!(import_packs %in% installed.packages()[,"Package"])]
if (length(notinstall) > 0) install.packages(notinstall)

sapply(load_packs, library, character.only = T)

wrap_dt <- function(dt, ...){
   para <- as.list(formals(DT::datatable))

   para$class <- "display compact"
   para$rownames = FALSE
   para$filter = "bottom"
   para$autoHideNavigation = FALSE
   para$escape = TRUE
   para$options <- list(scrollX = TRUE, scrollY = FALSE, autowidth = TRUE)
   input_list <- list(...)

   for (i in setdiff(names(input_list), "options")) {
      para[i] <- input_list[[i]]
   }

   for (i in names(input_list$options)) {
      para$options[i] <- input_list$options[[i]]
   }

   para$data <- dt
   do.call(what = datatable, para)
}
