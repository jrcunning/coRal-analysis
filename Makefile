all: JTB_ms/JTB_ms.pdf JTB_supp/JTB_supp.html

JTB_ms/JTB_ms.pdf: JTB_ms/JTB_ms.Rmd img/Fig1.png img/Fig2.png img/Fig3.png img/Fig4.png img/Fig5.png img/Fig6.png img/Fig7.png img/Fig8.png JTB_ms/library.bib
	R -e 'if(Sys.info()[["sysname"]]=="Darwin") { Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc") } else { Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio/bin/pandoc") }; rmarkdown::render("JTB_ms/JTB_ms.Rmd")'

img/Fig2.png: R/Fig2.R R/plot_steady_states.R R/run_steady_states.R R/run_coral_ss.R R/init_env.R R/def_pars.R
	R --vanilla < R/Fig2.R
	
img/Fig3.png: R/Fig3.R R/sensitivity.R R/run_coral_ss.R R/def_pars.R
	R --vanilla < R/Fig3.R
	
img/Fig4.png: R/Fig4.R R/run_coral.R R/run_coral_ss.R R/init_env.R R/def_pars.R
	R --vanilla < R/Fig4.R
	
img/Fig5.png: R/Fig5.R R/run_coral.R R/run_coral_ss.R R/init_env.R R/def_pars.R
	R --vanilla < R/Fig5.R
	
img/Fig6.png: R/Fig6.R R/run_coral.R R/run_coral_ss.R R/init_env.R R/def_pars.R
	R --vanilla < R/Fig6.R
	
img/Fig7.png: R/Fig7.R R/run_coral.R R/run_coral_ss.R R/init_env.R R/def_pars.R
	R --vanilla < R/Fig7.R
	
img/Fig8.png: R/Fig8.R R/run_coral.R R/run_coral_ss.R R/init_env.R R/def_pars.R
	R --vanilla < R/Fig8.R
	
JTB_supp/JTB_supp.html: JTB_supp/JTB_supp.Rmd JTB_supp/library.bib
	R -e 'if(Sys.info()[["sysname"]]=="Darwin") { Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc") } else { Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio/bin/pandoc") }; rmarkdown::render("JTB_supp/JTB_supp.Rmd")'
