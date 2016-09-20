JTB_ms/JTB_ms.pdf: JTB_ms/JTB_ms.Rmd img/Fig2.png img/model.png
	R -e 'if(Sys.info()[["sysname"]]=="Darwin") { Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc") } else { Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio/bin/pandoc") }; rmarkdown::render("JTB_ms/JTB_ms.Rmd")'

img/Fig2.png: R/Fig2_plot.R output/ss.Rdata
	R --vanilla < R/Fig2_plot.R

output/ss.Rdata: R/Fig2_steadystates.R R/run_coral.R R/init_env.R R/def_pars.R
	R --vanilla < R/Fig2_steadystates.R