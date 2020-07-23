
mainDir = getDirectory("Choose a main directory "); 
mainList = getFileList(mainDir); 
File.makeDirectory(mainDir+"Masks");
saveDir=mainDir+"Masks/"
Dialog.create("Create Masks on a Directory of Sub Directories");
Dialog.addString("image Suffix:", "ome");


Dialog.show();
imageSuffix = Dialog.getString() + ".tiff";


			
           		for (m=0; m<mainList.length; m++) { //clunky, loops thru all items in folder looking for image
		    		if (endsWith(mainList[m], imageSuffix)) { 
		   	 			
		   	 			open(mainDir+mainList[m]);
		   	 			im_title = getTitle();
			    		index = lastIndexOf(im_title, imageSuffix);
         				name = substring(im_title, 0, index);
						selectWindow(im_title);
						
						run("Gaussian Blur...", "sigma=3 stack");
						setAutoThreshold("Huang dark");
						setOption("BlackBackground", true);
						run("Convert to Mask", "method=Huang background=Dark black");
						run("Fill Holes", "stack");
						waitForUser;
						saveAs("tiff", saveDir+name+"_cell");

         				run("Set Measurements...", "area centroid perimeter shape stack redirect=None decimal=3");
          				run("Analyze Particles...", "size=1200-20000 display add stack");

       				 	array1 = newArray("0"); 
						for (i=1;i<roiManager("count");i++){ 
       				 	array1 = Array.concat(array1,i); 
        				Array.print(array1); 
						} 
						roiManager("select", array1); 
          
       					roiManager("save selected",  saveDir+name+"_roimask.zip");

						saveAs("Results", saveDir+name+"Tracker_Results.csv");
						run("Close");
						run("Close");
	

							    				
		    			}
           			}          						    

run("Close All");


