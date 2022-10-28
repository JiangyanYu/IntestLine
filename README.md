<img src="IntestLine-Logo.png" align="right" width=120 height=139 alt="" />

# IntestLine 
IntestLine is developed to map images of intestinal tissues prepared by Swiss-roll technique onto a line. It is suitable to process CODEX images or 10x spatial transcriptomics images, basically any images providing xy-coordinates.

# User manual
Tutorial on how to use the app can be found on Youtube: [https://www.youtube.com/watch?v=-QMrW_MKLPo](https://www.youtube.com/watch?v=-QMrW_MKLPo)

[Option 1: Use IntestLine via docker](#option1-docker)

[Option 2: Use IntestLine application implemented in FASTGenomics](#option2-fastgenomics)
A quick way to try IntestLine. Be aware that the configureation of the server is very limited (2 vCPUs @2.1GHz, 16 GB RAM). Therefore it will take a while to unroll a large dataset. For example, 26s for the demo dataset containing 10k cells, and xx min for a dataset with 150k cells.

[Option 3: Use IntestLine locally (Advanced users)](#option3-local-machine)

# Option1: Docker
1. Pull the docker image to your system
  ```
  docker pull altayyuzeir/intestline
  ```
2. Start the application
  ```
  docker run --name IntestLine -it -p 3838:3838 altayyuzeir/intestline bash
  ```
3. Once getting into the terminal, open the shiny app
  ```
  /usr/bin/run-shiny.sh shiny
  ```
4. Go to your web broswer and visit the port: http://localhost:3838/proxy/shiny
5. Do your analysis according to the steps described in option 2.
6. Ctrl+C in command terminal to close shiny app. exit - command to close the application

# Option2: FASTGenomics
Step 1: Upload data (Upper left panel in grey)\
1.1 Go to website https://beta.fastgenomics.org/a/intestline To use the application, unfortunately you need to login with an account. You can register a your own account, but to test the application you can also use our IntestLine account (Username: intestline@gmail.com Password: intestline).\
1.2 Upload CODEX-exported .csv file, containing x and y coordinates.\
1.3 If available, you can upload previously created backbone file.\
1.4 Select Region of Interest (ROI) by first choosing a center for the image.\
1.5 Then choose parameters a and b to highlight the elliptical area of interest.

Step 2: Select backbone points (Upper right panel in blue)\
2.1 Create a backbone by pressing with the mouse on the edge of the structures.\
2.2 **You must start from the center and go outward !**\
2.3 You can download the backbone selection for later use with the same ROI of the same image.\
2.4 You can preview the quality of the backbone selection. High quality backbone has a 3D spiral shape.

Step 3: Unroll (Bottom panel in green)\
3.1 Misprojected points can be filtered out by manipulating the angleCBS filter.\
3.2 Noisy points or contaminants can be eliminated by using the Z-score filtering.\
3.3 In both cases you can examine how many cells and from which regions are lost.

Step 4: Overlay parameters on the linear structure\
4.1 You can overlay expressed markers over the stretched image.\
4.2 This can be used to examine marker expression along the length and thickness of the organ.

# Option3: Local machine
1. Following libraries are needed to run the shiny app: igraph,tidyr,dplyr,ggplot2,plotly,ggtext,data.table,readr.
2. Download/clone the package from the IntestLine repository.
3. Run the app in docker/APP/app.R
4. Following the steps described in option 2.

# Cite us
Yuzeir, A.; Bejaran, D.; Grein, S.; Hasenauer, J.; Schlitzer, A.; Yu, J. IntestLine: A Shiny-Based Application to Map the Rolled Intestinal Tissue onto a Line; preprint; Bioinformatics, 2022. https://doi.org/10.1101/2022.10.26.513827.



