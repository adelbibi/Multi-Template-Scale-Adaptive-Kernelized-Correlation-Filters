"Multi-Template Scale-Adaptive Kernelized Correlation Filters" ICCVW2015

Authors: Adel Bibi and Bernard Ghanem.

Visit our group's website:
https://ivul.kaust.edu.sa/Pages/Home.aspx

Adel Bibi's website:
www.adelbibi.com

Email: adel.bibi@kaust.edu.sa
       bibiadel93@gmail.com


Bernard Ghanem's website:
http://www.bernardghanem.com/

Email: Bernard.Ghanem@kaust.edu.sa 

**************************************************
This MATLAB code implements a simple tracking pipeline based on the multi template scale
adpative kernelized correlation fitler (KCF_MTS).

It is free for research use. If you find it useful, please acknowledge the paper
above with a reference.
**************************************************


To run over 1 single video:

1- Open ('run_tracker.m').
2- Change line line 49 video = 'choose'.
3- Dump any video of OTB100/OTB50 into the the directory (videos).
4- Run the tracker.

The annotation files and the attributes for all the videos of OTB100 are avilable in
the directory (anno).


To run over the complete dataset:

1- Open ('run_tracker.m').
2- Change line line 49 video = 'all'.
3- Dump all the OTB100 or OTB50 videos into the the directory (videos).
4- Run the tracker.
5- The detailed results will be stored in (results) directory.


The code is also completey integratable with the OTB100 and OTB50 evaulation benchmarks.
To do so:

1- Move the complete traker directory to the (Trackers directory in the OTB evulation code).
The function is called (configTrackers.m) in the OTB evaulation code. It can be found here:
[1] http://cvlab.hanyang.ac.kr/tracker_benchmark/datasets.html
[2] https://sites.google.com/site/trackerbenchmark/benchmarks/v10

2- Add the following line to the list of trackers to be evualted over:
struct('name','KCF_MTSA','namePaper','KCF_MTSA')
Note: The code that will be run is through the evaulation is (run_KCF_MTSA.m).
It's the same exact code with the same parameters but has been changed to the standard format.


**************************************************

Referrences:
[1] Henriques, João F., et al. "High-speed tracking with kernelized correlation filters." Pattern Analysis and Machine Intelligence, IEEE Transactions on 37.3 (2015): 583-596.
[2] Henriques, Joao F., et al. "Exploiting the circulant structure of tracking-by-detection with kernels." Computer Vision–ECCV 2012. Springer Berlin Heidelberg, 2012. 702-715.
[3] Wu, Yi, Jongwoo Lim, and Ming-Hsuan Yang. "Online object tracking: A benchmark." Proceedings of the IEEE conference on computer vision and pattern recognition. 2013.

A complete list of references can be found in the paper, which can be found here

https://ivul.kaust.edu.sa/Pages/Pub-Adaptive-Kernelized-Correlation-Filters.aspx
http://www.adelbibi.com/papers/ICCVW2015/paper.pdf


