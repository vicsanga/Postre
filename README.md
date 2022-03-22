# Postre
Software developed to predict and annotate the impact of Structural Variants affecting both the coding and non-coding genome.
<ul>
      <li><a href="#ExplanationPostre">What is Postre?</a></li>
      <li><a href="#UsingPostre">How to use Postre?</a></li>
      <li><a href="#Installation">How to install and run Postre</a></li>
</ul>
<h2 id="ExplanationPostre"> <b>What is Postre?</b> </h2>

 <div>
Postre is a software developed to <b style='color:#1D3354;'>predict the pathogenic impact of Structural Variants (SVs)</b>. Postre aims to identify the genes responsible for the disease together with the pathological mechanism. In parallel, it determines the SV associated pathogenic events <b style='color:#1D3354;'>in space and time</b>. Postre is able to handle <b style='color:#1D3354;'>long-range (enhancer driven) and coding pathogenic events</b>.
 <br> <br>
</div>

Check the infographic displayed below to see a representation of Postre functionality. <b>Click the image to expand it and see it in high resolution!</b>

![Postre Diagram](https://github.com/vicsanga/Postre/blob/main/Postre_app/www/infografia_test_600ppp.png?raw=true)

<h2 id="UsingPostre">How to use Postre?</h2>

A quick tutorial showing main Postre features and explaining its usage is provided in the video below. 

& Video youtube

<br><br>
<h2 id="Installation">How to install and run Postre</h2>

This is the explanation to use Postre locally, using your device resources, without the need to depend on the web option availability. At the same time, you will also be running Postre latest version.

<b>Eventhough you may not have computational skills, do not be afraid!. Postre installation is very easy.</b> On top of that, once Postre is installed, running and using it is as simple as any other desktop application. Thanks to its user friendly graphical interface.

Postre is built with <a href="https://shiny.rstudio.com/" target="_blank">Shiny</a> framework.
Thus, to run Postre you only require R (version >=3.5.0) and a couple of R libraries.

<h3>1 Installing R </h3>
There is plenty of available information in the internet to do this depending on the OS (Windows, Mac, Linux etc.). If you need help, different tutorials are provided here: 
<br><br>
<ul>
<li><a href="https://www.youtube.com/watch?v=NZxSA80lF1I" target="_blank">Windows Youtube Tutorial </a></li>
<li><a href="https://www.youtube.com/watch?v=LanBozXJjOk" target="_blank">Mac Youtube Tutorial </a></li>
<li><a href="https://www.youtube.com/watch?v=iN0UZ43G6GE"target="_blank">Ubuntu Youtube Tutorial </a></li>
<li><a href="https://www.earthdatascience.org/courses/earth-analytics/document-your-science/setup-r-rstudio/">Web Tutorial for installing R and R-Studio in Windows, Mac or Linux <a/></li>
</ul>

<h3>2. Running Postre Latest Version</h3>     
Once R  is installed, just execute the following 2 instructions in a R console to run Postre latest version:
<br><br>

```R
##Instruction to manage dependencies and run postre
source("https://raw.githubusercontent.com/vicsanga/Postre/main/Postre_wrapper.R")
```

The above comand installs and loads all the required libraries for Postre. If some libraries are missing R may ask for permission to install them, confirm it. The first time that you run Postre this action will probably take more time.

If you are not sure about how to do it. You can find me initializing Postre in the video below!
<br><br>

YOUTUBE VIDEO INITIALIZING POSTRE
      
<br><br><br><br>
