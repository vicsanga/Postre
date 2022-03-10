# Postre
Software developed to predict and annotate the impact of Structural Variants affecting both the coding and non-coding genome.
<ul>
      <li><a href="#ExplanationPostre">What is Postre?</a></li>
      <li><a href="#UsingPostre">Using Postre</a></li>
      <ul>
         <li><a href="#Installation">Local usage: How to install and run Postre on your computer</a></li>
         <li><a href="#OnlinePostre">Online usage: Running postre on the web</a></li>
      </ul>
</ul>
<h2 id="ExplanationPostre"> <b>What is Postre?</b> </h2>

 <div>
Postre is a software developed to <b style='color:#1D3354;'>predict the pathogenic impact of Structural Variants (SVs)</b>. Postre aims to identify the genes responsible for the disease together with the pathological mechanism. In parallel, it determines the SV associated pathogenic events <b style='color:#1D3354;'>in space and time</b>. Postre is able to handle <b style='color:#1D3354;'>long-range (enhancer driven) and direct pathogenic events</b>.
 <br> <br>
</div>

![Postre Diagram](https://github.com/vicsanga/Postre/blob/main/Postre_app/www/WhatIsPostre.png?raw=true)

& Video youtube


<h2 id="UsingPostre">Using Postre</h2>

Postre can be used in two different ways:

 <ul>
   <li><a href="#Installation">Local usage: How to install and run Postre on your computer</a></li>
   <li><a href="#OnlinePostre">Online usage: Running postre on the web</a></li>
 </ul>
      
<i>Online option may not ve available if users usage has already exceeded the cloud suscribed services (e.g. maximum hours usage).</i>      

<h2 id="Installation">Local usage: How to install and run Postre on your computer</h2>

This is the explanation to use Postre locally, using your device resources, without the need to depend on the web option availability. At the same time, you will be running Postre latest version. Before going into detail, <b>eventhough you may not have computational skills, do not be afraid!. Postre installation is very easy.</b> On top of that, once you (or the IT colleague) install Postre, running and using it is as easy as any other desktop application. Thanks to its user friendly graphical interface.

Postre is built with <a href="https://shiny.rstudio.com/" target="_blank">Shiny</a> framework.
Thus, to run Postre you only require R & R-Studio, and a couple of R libraries.

<h3>1. Installing R and R-Studio </h3>
There is plenty of available information in the internet to do this depending on the OS (Windows, Mac, Linux etc.). If you need help, different tutorials are provided here: 
<ul>
<li><a href="https://www.youtube.com/watch?v=NZxSA80lF1I" target="_blank">Windows Youtube Tutorial </a></li>
<li><a href="https://www.youtube.com/watch?v=LanBozXJjOk" target="_blank">Mac Youtube Tutorial </a></li>
<li><a href="https://www.youtube.com/watch?v=iN0UZ43G6GE"target="_blank">Ubuntu Youtube Tutorial </a></li>
<li><a href="https://www.earthdatascience.org/courses/earth-analytics/document-your-science/setup-r-rstudio/">Web Tutorial for installing R and R-Studio in Windows, Mac or Linux <a/></li>
</ul>

<h3>2. Running Postre</h3>      
It is very straightforward to run Postre latest version.  It is only required to execute the instructions provided below in R. First, open R studio, afterwards run the following 2 instructions: 
<br><br>

```R
##Instruction 1: Managing & Loading Dependencies
source("https://raw.githubusercontent.com/vicsanga/Postre/main/Managing_Postre_Dependencies.R")

##Instruction 2: Running Postre
runGitHub("Postre", "vicsanga", subdir = "Postre_app/", launch.browser = TRUE)
```

If you want to know more about the different instructions:

<ol>
<li><b>Instruction 1 </b>: Installs and loads all the required libraries for Postre. If some libraries are missing R may ask for permission to install them, confirm it.The first time that you run Postre this action will probably take more time. Instruction 1 must always be executed, because if a new library is added to Postre it will be considered in this first step.</li>

<li><b>Instruction 2 </b>: Runs Postre latest version in your web browser.</li>

</ol>
You can just copy paste and run those instructions at once in R. You can find me doing that in the video below!

<br><br>

YOUTUBE VIDEO INITIALIZING POSTRE
      

<br><br>

<h2 id="OnlinePostre">Online usage: Running postre on the web</h2>
Postre latest version is also hosted in the cloud to offer an online usage. However, given limitations of the currently hired cloud service, <b>it may not be available if the cloud monthly suscribed services have been consumed, (e.g. maximum usage hours)</b>. If you find Postre unavailable online, go for <a href="#Installation">Local usage: How to install and run Postre in your computer</a>. 

<br><br>
You can find Postre online here: <a href="https://svradalab.shinyapps.io/postre_app/">Postre online service</a>.

<br><br><br><br>
