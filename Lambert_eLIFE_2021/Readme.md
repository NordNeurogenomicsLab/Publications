### Manuscript

[Jason T. Lambert, Linda Su-Feher, Karol Cichewicz, Tracy L. Warren, Iva Zdilar, Yurong Wang, Kenneth J. Lim, Jessica Haigh, Sarah J. Morse, Cesar P. Canales, Tyler W. Stradleigh, Erika Castillo, Viktoria Haghani, Spencer Moss, Hannah Parolini, Diana Quintero, Diwash Shrestha, Daniel Vogt, Leah C. Byrne, Alex S. Nord (2021).
 **Parallel functional testing identifies enhancers active in early postnatal mouse brain.**](https://www.biorxiv.org/content/10.1101/2021.01.15.426772v3)


### R Markdown html analysis reports referenced in the manuscript

1. [MPRA](https://nordneurogenomicslab.github.io/STAR408/)     
2. [miniMPRA](https://nordneurogenomicslab.github.io/miniMPRA/)

### Links to the analysis repository
1. [MPRA repository](https://github.com/NordNeurogenomicsLab/STAR408)
2. [miniMPRA repository](https://github.com/NordNeurogenomicsLab/miniMPRA)   

### Links to the analysis Docker image repositories
1. [The major MPRA analysis image](https://hub.docker.com/repository/docker/kcbio/lambert_elife_2021_star408)
2. [miniMPRA analysis image](https://hub.docker.com/repository/docker/kcbio/mini_mpra)

Docker images are based on the rocker/tidyverse image, running R 4.1.0. All necessary packages are preinstalled. Docker images launch RStudio in a web browser.

#### Launch Docker images from a terminal:
1. docker run -d --rm --name STAR408 -p 8787:8787 kcbio/lambert_elife_2021_star408:final_Rmd
2. docker run -d --rm --name miniMPRA -p 8787:8787 kcbio/mini_mpra:final_image

Go to http://localhost:8787/ in your webbrowser
Specify -p 8788:8787 if you want to run both images concurrently. 
