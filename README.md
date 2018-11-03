# AG3line:Active Grouping and Geometry-Gradient Combined Validation for Line Segment Extraction

## Introduction
This is the implementation of the proposed line segment detector named AG3line: "Active Grouping and Geometry-Gradient Combined Validation for Line Segment Extraction". The extraction result will be shown in figure after running the code in 'AG3lineDemo.m', and the only parameter in the code is the name of picture to be used in folder "pic". The folder "pic" contains 13 images for experiment, and about 100 images can be found on the [website](http://www.elderlab.yorku.ca/resources/york-urban-line-segment-database-information/), the York Urban database. Also, the experimental results of AG3line, [LSD](http://www.ipol.im/pub/art/2012/gjmr-lsd/), [Linelet](https://github.com/NamgyuCho/Linelet-code-and-YorkUrban-LineSegment-DB) and [EDLines](http://ceng.anadolu.edu.tr/cv/EDLines/) on the York Urban database are saved in folder "YorkUrbanRes", and the comparing of the extraction result can be plotted by the code in "comparedraw.m".  
Since this is the experimental vision coded by MATLAB, the running speed remains to be improved. It takes about 10-20 seconds to deal with the image in size of 640 x 480, and the image over 1000 in width or height is not recommended in experiment since it may take a long time. Good luck!
## Some results
Some experimental results are listed below. In addition, the results of other 100 images in York Urban databse are saved in "YorkUrbanRes" and they can be ploted by the code in "YorkUrbanRes\seeres.m".  

|||||
|----|----|----|----|
|image| PPHT | LSD | EDLines |
|![](pic/testimg/010.png "005") |![](pic/result/tower1-ppht.png "005")|![](pic/result/tower1-LSD.jpg "005")|![](pic/result/tower1-EDLines.jpg "005")|
|| Linelet | MCMLSD | AG3line |
||![](pic/result/tower1-Linelet.jpg "005")|![](pic/result/tower1-mcm.jpg "005")|![](pic/result/tower1-AG3line.jpg "005")|
|image| PPHT | LSD | EDLines |
|![](pic/testimg/tower11.jpg "005")|![](pic/result/tower2-ppht.png "005")|![](pic/result/tower2-LSD.jpg "005")|![](pic/result/tower2-EDLines.jpg "005")|
|| Linelet | MCMLSD | AG3line |
||![](pic/result/tower2-Linelet.jpg "005")|![](pic/result/tower2-mcm.jpg "005")|![](pic/result/tower2-AG3line.jpg "005")|
|image| PPHT | LSD | EDLines |
|![](pic/testimg/P1020856.jpg "005")|![](pic/result/facade-ppht.jpg "005")|![](pic/result/facade-LSD.jpg "005")|![](pic/result/facade-EDLines.jpg "005")|
|| Linelet | MCMLSD | AG3line |
||![](pic/result/facade-Linelet.jpg "005")|![](pic/result/facade-mcm.jpg "005")|![](pic/result/facade-AG3line.jpg "005")|
|image| PPHT | LSD | EDLines |
|![](pic/testimg/P1080091.jpg "005")|![](pic/result/york3-ppht.jpg "005")|![](pic/result/your3-LSD.jpg "005")|![](pic/result/york3-EDLines.jpg "005")|
|| Linelet | MCMLSD | AG3line |
||![](pic/result/york3-Linelet.jpg "005")|![](pic/result/york3-mcm.jpg "005")|![](pic/result/york3-AG3line.jpg "005")|
|image| PPHT | LSD | EDLines |
|![](pic/testimg/P1080005.jpg "005")|![](pic/result/york1-ppht.jpg "005")|![](pic/result/york1-LSD.jpg "005")|![](pic/result/york1-EDLines.jpg "005")|
|| Linelet | MCMLSD | AG3line |
||![](pic/result/york1-Linelet.jpg "005")|![](pic/result/york1-mcm.jpg "005")|![](pic/result/york1-AG3line.jpg "005")|







