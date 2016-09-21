Rmd/JTB_ms.pdf: Rmd/JTB_ms.Rmd img/Fig2.png img/model.png
	R -e 'Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc");rmarkdown::render("Rmd/JTB_ms.Rmd")'

img/Fig2.png: R/Fig2.R
	R --vanilla < R/Fig2.R