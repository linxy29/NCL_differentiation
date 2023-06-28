knitr::opts_chunk$set(echo = TRUE, results = "markup", error = TRUE, fill = TRUE, out.width = "100%")
options(knitr.kable.NA = "", width = 1000)

load_packs  <- c("knitr","data.table", "DT", "stringr", "ggplot2")
require_packs <- c("tools", "vioplot", "ggdag", "htmlTable" )
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

autolayout <- function(n, x, min = "row", byrow = F, nr, nc){
  if (missing(n) ) n <- length(x)

  if (missing(nr) & missing(nc)) {
    n1 <- ceiling(sqrt(n))
    n2 <- ceiling(n / n1)

    if (min == "row") {
      nr <- min(c(n1, n2))
      nc <- max(c(n1, n2))
    } else {
      nc <- min(c(n1, n2))
      nr <- max(c(n1, n2))
    }
  } else {
    if (missing(nr)) {
      nr <- ceiling(n / nc)
    } else {
      nc <- ceiling(n / nr)
    }
  }
  layout(matrix(1:(nc * nr), byrow = byrow, ncol = nc, nrow = nr))
}

thepar <- list(bty = "l",  mar = c(3,3,.5,.5),mgp = c(1.25,.25,0), oma = c(0,0,2,0), font.lab = 2, tck = -0.02, cex = 1)

subset_raw <- function(dt = dt_raw, cible = "ct"){
  if (!cible %in% c("ct", "delta_ct", "delta_delta_ct", "fge", "log2_fge", "relative_count")) stop("cible must be ct, delta_ct or delta_delta_ct or relative_count")
  if (cible == "ct") return(dt)
  if (cible == "delta_ct") {
    dt_sub <- dt[as.character(gene) != as.character(ct_control_gene)]
    dt_sub[, gene := droplevels(gene)]
    return(dt_sub)
  }

  if (cible %in% c( "delta_delta_ct", "fge", "log2_fge")) {
    dt_sub <- dt[as.character(gene) != as.character(ct_control_gene) & day != 0,
                 list(delta_delta_ct = mean(delta_delta_ct, na.rm  = T),
                      fge = mean(fge, na.rm = T),
                      log2_fge = mean(log2_fge, na.rm = T)), by = .(gene, sb, sb_conc, sb_day, exp, clone, day, CHIR, NOTO_mRNA, extra_replicate)]
    dt_sub[, gene := droplevels(gene)]
    return(dt_sub)
  }
  if (cible == "relative_count") {
    return(dt)
  }
}

extract_hyperlink_20220506_20220707 <- function(current){
  output <- list.files(file.path("..", "reports"))
  output <- output[stringr::str_detect(output, "ipspine_SB_treatment_20220506_20220707")]
  output <- output[stringr::str_detect(output, current, negate = T)]
  links <- paste0("[", output, "](",output,")")
  texte <- paste0("Other reports about the analyse of SB treatement data of 20220506 and 20220707 follow: ", knitr::combine_words(links), ".")

  return(texte)
}

extract_hyperlink_20220811 <- function(current){
  output <- list.files(file.path("..", "reports"))
  output <- output[stringr::str_detect(output, "ipspine_SB_treatment_20220811")]
  output <- output[stringr::str_detect(output, current, negate = T)]
  links <- paste0("[", output, "](",output,")")
  texte <- paste0("Other reports about the analyse of SB treatement data of 20220811 follow: ", knitr::combine_words(links), ".")

  return(texte)
}

extract_hyperlink_20220712 <- function(current){
  output <- list.files(file.path("..", "reports"))
  output <- output[stringr::str_detect(output, "ipspine_noSB_treatment_20220712")]
  output <- output[stringr::str_detect(output, current, negate = T)]
  links <- paste0("[", output, "](",output,")")
  texte <- paste0("Other reports about the analyse of no SB treatement data of 20220712 follow: ", knitr::combine_words(links), ".")

  return(texte)
}

extract_hyperlink_20220817 <- function(current){
  output <- list.files(file.path("..", "reports"))
  output <- output[stringr::str_detect(output, "ipspine_quantification_20220817")]
  output <- output[stringr::str_detect(output, current, negate = T)]
  links <- paste0("[", output, "](",output,")")
  texte <- paste0("Other reports about the analyse of quantification data of 20220817 follow: ", knitr::combine_words(links), ".")

  return(texte)
}

extract_hyperlink_20230418 <- function(current){
  output <- list.files(file.path("..", "reports"))
  output <- output[stringr::str_detect(output, "ipspine_noSB_treatment_20230418")]
  output <- output[stringr::str_detect(output, current, negate = T)]
  links <- paste0("[", output, "](",output,")")
  texte <- paste0("Other reports about the analyse of no SB treatement data of 20230418 follow: ", knitr::combine_words(links), ".")

  return(texte)
}

col_quali2 <- c("#33a02c", "#1f78b4")
col_quali3 <- c("#1b9e77", "#d95f02", "#7570b3")
col_quali4 <- c("#a6cee3", "#1f78b4", "#b2df8a","#33a02c")
col_quali9 <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf","#999999")

col_seq9 <- c("#fff7ec","#fee8c8","#fdd49e","#fdbb84","#fc8d59","#ef6548","#d7301f","#b30000","#7f0000")

# Fonctiuns for run reports from external script

preambule <- function(titre){
  texte <- paste0(c(
    "---\n",
    "title: ", titre, "\n",
    "author: Antoine Humeau","\n",
    "date: ", as.character(Sys.time()), "\n",
    "output: ", "\n",
    "  html_document:", "\n",
    "    code_folding: hide", "\n",
    "    fig_caption: yes", "\n",
    "    highlight: tango", "\n",
    "    number_section: yes", "\n",
    "    theme: cosmo", "\n",
    "    toc: yes", "\n",
    "    toc_float: yes", "\n",
    "    toc_depth: 4", "\n",
    "    keep_md: no", "\n",
    "    df_print: paged", "\n",
    "---", "\n",
    "<style>", "\n",
    ".main-container {width: 80%;max-width: inherit;}", "\n",
    "</style>", "\n"

  ), collapse = "")

  readLines(textConnection(texte))
}

generate_preambule_rmd <- function(titre, rmd_file = "temporaire.Rmd", erase = F) {
  book_header <- preambule(titre = titre)

  if (length(grep(rmd_file, list.files())) > 0) {
    if (erase) {
      warning(paste(rmd_file, "already exists and was recreated"))
    } else {
      stop(paste(rmd_file, "already exists. You can paste erase =TRUE if you want to replace the old document"))
    }
  }
  write(book_header, file = rmd_file )
}

add_template_to_rmd <- function(template, rmd_file = "temporaire.Rmd", titre, titre_level, parameters = NULL){
  texte <- NULL
  if (!missing(titre)) {
    if (missing(titre_level)) titre_level = 1
    texte <- paste(paste0(rep("#", titre_level), collapse = "" ), titre)
  }

  if (!is.null(parameters)) {
    for (k in names(parameters)) {
      tmp <- parameters[[k]]
      new_object <- paste0("the", k)
      assign(x = new_object, value = tmp, envir = .GlobalEnv)
    }
  }
  texte <- c(texte, readLines(template), sep = "\n")

  write(texte, sep = "\n", file = rmd_file, append = T)
}

generate_new_rmd <- function(titre, rmd_file = "temporaire.Rmd", erase = FALSE,
                             templates_add, templates_dir = "templates"){

  generate_preambule_rmd(titre = titre, rmd_file = rmd_file, erase = erase)

  for (i in 1:length(templates_add)) {
    print(paste("Run the", names(templates_add)[i], "template"))

    theargs <- templates_add[[i]]
    theargs$template <- names(templates_add)[[i]]

    theargs_add_template <- formals(add_template_to_rmd)
    matchs_args <- match(names(theargs_add_template), names(theargs))

    theargs_add_template[!is.na(matchs_args)] <- theargs[na.omit(matchs_args)]
    theargs_add_template$rmd_file <- rmd_file

    do.call(add_template_to_rmd, c(theargs_add_template))
  }
}

generate_new_report <- function(titre, rmd_file = "temporaire.Rmd", templates_add, templates_dir = "templates",
                                erase = FALSE, clean = T, clean_rmd = F, quiet = F){

  theargs_generate_new_rmd <- formals(generate_new_rmd)

  input_params <- as.list(match.call(generate_new_report))

  matchs_args <- match(names(theargs_generate_new_rmd), names(input_params))
  theargs_generate_new_rmd[!is.na(matchs_args)] <- input_params[na.omit(matchs_args)]

  do.call(generate_new_rmd, c(theargs_generate_new_rmd))

  rmarkdown::render(rmd_file, encoding = "UTF-8", output_dir = "../reports", envir = new.env(), clean = clean, quiet = quiet)

  if (clean_rmd) {
    unlink(rmd_file)
    warning(paste(rmd_file, "was deleted after the html creation"))
  }
}

select_template <- function(x, templates_dir = "templates", templates_prefix = ""){
  if (templates_prefix != "") x <- paste0(templates_prefix, x)
  x[tools::file_ext(x) == ""] <- paste0(x[tools::file_ext(x) == ""], ".Rmd")
  if (templates_dir != "") x <- file.path(templates_dir, x)
  return(x)
}

wrap_generate_new_report_analyses <- function(parameters){
  names(parameters$templates_add) <- select_template(names(parameters$templates_add), templates_dir = parameters$templates_dir, templates_prefix =  parameters$templates_prefix)
  theargs <- formals(generate_new_report)
  matchs_args <- match(names(theargs), names(parameters))
  theargs[!is.na(matchs_args)] <- parameters[na.omit(matchs_args)]

  do.call(generate_new_report, c(theargs))
}
