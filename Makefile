all: ms/ms.pdf supp/supp.html

ms/ms.pdf: ms/ms.Rmd ms/library.bib img/*.png
	R -e 'if(Sys.info()[["sysname"]]=="Darwin") { Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc") } else { Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio/bin/pandoc") }; rmarkdown::render("ms/ms.Rmd")'

img/Fig2.png: R/Fig2.R 
	R --vanilla < R/Fig2.R
	
img/Fig3.png: R/Fig3.R R/plot_steady_states.R R/run_steady_states.R
	R --vanilla < R/Fig3.R
	
img/Fig4.png: R/Fig4.R R/sensitivity.R
	R --vanilla < R/Fig4.R
	
img/Fig5.png: R/Fig5.R
	R --vanilla < R/Fig5.R
	
img/Fig6.png: R/Fig6.R
	R --vanilla < R/Fig6.R
	
img/Fig7.png: R/Fig7.R
	R --vanilla < R/Fig7.R
	
img/Fig8.png: R/Fig8.R
	R --vanilla < R/Fig8.R
	
img/Fig9.png: R/Fig9.R
	R --vanilla < R/Fig9.R
	
supp/supp.html: supp/supp.Rmd supp/library.bib
	R -e 'if(Sys.info()[["sysname"]]=="Darwin") { Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc") } else { Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio/bin/pandoc") }; rmarkdown::render("supp/supp.Rmd")'
