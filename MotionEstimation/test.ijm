//open stack (Kennt sich einer aus wie man den Pfad nicht absolut setzen muss?)
run("AVI...", "open=[D:\\Dropbox\\Eclipse_Workspace\\ABV_Bewegungsanalyse\\MotionEstimation\\b_testmaterial2.avi] first=0 last=92 use");
run("32-bit");

//motion estimation
run("MotionVectorField ");