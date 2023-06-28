source(file.path("..", "R_functions", "analyse_prelude.R"))

# No sb treatement  ----

### 20230610 ----

params_des_nosb <- list(titre = "RT-qPCR Figures 1B and supplementary",
                        erase = T, quiet = TRUE,
                        templates_add = list(intro = list(titre = "load", titre_level = 1, parameters = list(hyperlink_curent = "data_description")),
                                             vioplot_call_1outcome = list(titre = "Outcome distribution by experimental variables", titre_level = 1),
                                             figs = list(titre = "Fold-change at day x relative to day 0 by gene and clone for NOTO_mRNA", titre_level = 1)
                        ),
                        templates_dir = file.path("templates", "noSB_treatment"),
                        templates_prefix = "ipspine_noSB_treatment_20230610_analyses_template_", 
                        
                        rmd_file = "ipspine_noSB_treatment_20230610_data_description.Rmd")

wrap_generate_new_report_analyses(params_des_nosb)



