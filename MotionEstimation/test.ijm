//Lukas
open("D:\\Dropbox\\Eclipse_Workspace\\ABV_Bewegungsanalyse\\MotionEstimation\\TIFF_sequence\\img_stack.tif");

//Laura 
//run("AVI...", "open=[//Users//laura//Documents//GitHub//ABV_Bewegungsanalyse//MotionEstimation//b_testmaterial2.avi] first=0 last=92 use");
//run("AVI...", "open=[//Users//laura//Google Drive//Bachelorarbeit//Testmaterial//einPixel.avi] first=0 last=92 use");

run("32-bit");

//motion estimation
run("MotionVectorField ");