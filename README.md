<img src="IntestLine-Logo.png" align="right" width=120 height=139 alt="" />

# IntestLine 
IntestLine is developed to map images of intestinal tissues prepared by Swiss-roll technique onto a line. It is suitable to process CODEX images or 10x spatial transcriptomics images, basically any images providing xy-coordinates.

# User manual
Tutorial on how to use the app can be found on Youtube: [https://www.youtube.com/watch?v=DLZ-4taQk-s](https://www.youtube.com/watch?v=DLZ-4taQk-s)

Option 1: Use IntestLine via docker
1. Pull the docker image to your system

  ```
  docker pull altayyuzeir/intestline
  ```
  
3. Open command terminal
4. cd C:\Users\altay\Altay\Lectures_Biochemustry\5_Lab_rotation_2\ - set desired work folder with cd command
5. docker pull altayyuzeir/intestline - to pull the image
6. docker run --name IntestLine -it -p 3838:3838 altayyuzeir/lab-rotation-2:fg-v0.60 bash - start the application
7. /usr/bin/run-shiny.sh shiny - this command opens the shiny app
8. visit: http://localhost:3838/proxy/shiny on the default browser
9. do your analysis
10. Ctrl+C in command terminal to close shiny app
11. exit - command to close the application

[Option 2: Use IntestLine application implemented in FASTGenomics](#option2-fastgenomics)\
A quick way to try IntestLine. Be aware that the configureation of the server is very limited (). therefore it will take a while to unroll a large dataset.

Option 3: Use IntestLine locally

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

# Cite us

This package is aiming to unroll intestinal images from CODEX.\
This script is inspired by Parigi et al (https://www.nature.com/articles/s41467-022-28497-0). \
The method is implemented in a shiny-based web application: https://beta.fastgenomics.org/a/intestline


