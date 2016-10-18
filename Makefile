JTB_ms/JTB_ms.pdf: JTB_ms/JTB_ms.Rmd img/Fig1.png img/Fig2.png img/Fig3.png JTB_ms/library.bib
	R -e 'if(Sys.info()[["sysname"]]=="Darwin") { Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc") } else { Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio/bin/pandoc") }; rmarkdown::render("JTB_ms/JTB_ms.Rmd")'

img/Fig2.png: R/Fig2.R R/plot_steady_states.R R/run_steady_states.R R/run_coral_ss.R R/init_env.R R/def_pars.R
	R --vanilla < R/Fig2.R
	
img/Fig3.png: R/Fig3.R R/sensitivity.R R/run_coral_ss.R R/def_pars.R
	R --vanilla < R/Fig3.R