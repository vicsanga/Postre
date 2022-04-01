# POSTRE
Software developed to predict and annotate the impact of Structural Variants affecting both the coding and non-coding genome.
<ul>
      <li><a href="#ExplanationPOSTRE">What is POSTRE?</a></li>
      <li><a href="#UsingPOSTRE">How to use POSTRE?</a></li>
      <li><a href="#Installation">How to install and run POSTRE</a></li>
</ul>
<h2 id="ExplanationPOSTRE"> <b>What is POSTRE?</b> </h2>

 <div>
POSTRE is a software developed to <b style='color:#1D3354;'>predict the pathogenic impact of Structural Variants (SVs)</b>. POSTRE aims to identify the genes responsible for the disease together with the pathological mechanism. In parallel, it determines the SV associated pathogenic events <b style='color:#1D3354;'>in space and time</b>. POSTRE is able to handle <b style='color:#1D3354;'>long-range (enhancer driven) and coding pathogenic events</b>.
 <br> <br>
</div>

Check the infographic displayed below to see a representation of POSTRE functionality. <b>Click the image to expand it.</b>

![Postre Diagram](https://github.com/vicsanga/Postre/blob/main/Postre_app/www/infografia.png?raw=true)

<h2 id="UsingPOSTRE">How to use POSTRE?</h2>

A quick tutorial showing main POSTRE features and explaining its usage is provided here. Click on the image below to see the video in Youtube. Reproduce it in Full Screen and High Quality (1080p) for optimal visualitzation. 

[![POSTRE Tutorial](https://github.com/vicsanga/Postre/blob/main/Postre_app/www/ImagenParaGithub_Tutorial.png?raw=true)](https://www.youtube.com/watch?v=SeR6vD0wPrE "POSTRE Tutorial")

<br><br>
<h2 id="Installation">How to install and run POSTRE</h2>

This is the explanation to use POSTRE locally, using your device resources, without the need to depend on the web option availability. At the same time, you will also be running POSTRE latest version.

<b>Eventhough you may not have computational skills, do not be afraid!. POSTRE installation is very easy.</b> On top of that, once POSTRE is installed, running and using it is as simple as any other desktop application. Thanks to its user friendly graphical interface.

POSTRE is built with <a href="https://shiny.rstudio.com/" target="_blank">Shiny</a> framework.
Thus, to run POSTRE you only require R (version >=3.5.0) and a couple of R libraries.

<h3>1 Installing R </h3>

To run Postre R version >=3.5.0 is required.

There is plenty of available information in the internet about how to install R depending on the OS (Windows, Mac, Linux etc.). Usually upon R installation people also install R-Studio, which is an integrated development environment to write R code. As a result, most tutorials explain how to install both R and R-Studio. But, to run Postre R-Studio is not necessary.  If you need help for R download, different tutorials are provided here: 
<br><br>
<ul>
<li><a href="https://www.youtube.com/watch?v=NZxSA80lF1I" target="_blank">Windows Youtube Tutorial </a></li>
<li><a href="https://www.youtube.com/watch?v=LanBozXJjOk" target="_blank">Mac Youtube Tutorial </a></li>
<li><a href="https://www.youtube.com/watch?v=iN0UZ43G6GE"target="_blank">Ubuntu Youtube Tutorial </a></li>
<li><a href="https://www.earthdatascience.org/courses/earth-analytics/document-your-science/setup-r-rstudio/">Web Tutorial for installing R and R-Studio in Windows, Mac or Linux <a/></li>
</ul>

<h3>2. Running POSTRE Latest Version</h3>     
Once R  is installed, just execute the following instruction in a R console to run POSTRE latest version:
<br><br>

```R
source("https://raw.githubusercontent.com/vicsanga/Postre/main/Postre_wrapper.R")
```

The above instruction installs and loads all the required libraries for POSTRE. If some libraries are missing R may ask for permission to install them, confirm it. The first time that you run POSTRE this action will probably take more time.

If you are not sure about how to do it. You can find me initializing POSTRE from R in the video below!
<br><br>

[![POSTRE running](https://github.com/vicsanga/Postre/blob/main/Postre_app/www/ImageGithub_HowToRunPostre.png?raw=true)](https://www.youtube.com/watch?v=FCYitDfbnXk "POSTRE Running")
      
<br><br><br><br>
