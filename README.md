# Postre
Software developed to predict and annotate the impact of Structural Variants affecting both the coding and non-coding genome.
<ul>
      <li><a href="#ExplanationPostre">What is Postre?</a></li>
      <li><a href="#Installation">How to install and run Postre?</a></li>
</ul>
<h2 id="ExplanationPostre"> <b>What is Postre?</b> </h2>

 <div>
Postre is a software developed to <b style='color:#1D3354;'>predict the pathogenic impact of Structural Variants (SVs)</b>. Postre aims to identify the genes responsible for the disease together with the pathological mechanism. In parallel, it determines the SV associated pathogenic events <b style='color:#1D3354;'>in space and time</b>. Postre is able to handle <b style='color:#1D3354;'>long-range (enhancer driven) and direct pathogenic events</b>.
 <br> <br>
</div>

IMAGEN

& Video youtube



<h2 id="Installation">How to install and run Postre?</h2>

Postre is built with <a href="https://shiny.rstudio.com/">Shiny</a> framework.
Thus, to run Postre you only need to install R & R-Studio, and a couple of R libraries.

<h3>1. Installing R and R-Studio </h3>
There is plenty of available information in the internet to do this depending on the OS(Windows, Mac, Linux etc.). If you need help, different tutorials are provided here: <a href="https://www.earthdatascience.org/courses/earth-analytics/document-your-science/setup-r-rstudio/">Tutorial for installing R and R-Studio in Windows, Mac or Linux <a/>

<h3>2. Running Postre</h3>      
It is very straightforward to run Postre latest version. It is only required to execute the 3 instructions provided below in R. First, open R studio, afterwards run the following 3 instructions: 

```R
if (!require("pacman")){install.packages("pacman")} ##Instruction 1: Checking pacman installed or installing
source()##Instruction 2: Installing and loading Postre required libraries
##Instruction 3: Running Postre (latest version)
```
Instruction 2 and 3 are always mandatory. 
If you want to know more about the different instructions. <b>Instruction 1 </b> checks if a library manager called pacman is installed, if it is not, then pacman will be (if R asks you, confirm the installation of the package). The <b>Instruction 2 </b> loads and installs all the required libraries for Postre, taking profit of pacman functionality (the first time that you run Postre this action will take more time, since the different libraries will be installed). The <b>Instruction 3 </b> will run and open Postre in your web browser.


You can just copy paste and run those instructions at once in R. You can find me doing that in the video below!

YOUTUBE VIDEO INITIALIZING POSTRE
      
      
Running Postre on the web?

To run it

```
To install Postre
source(whtvr.R)
```

