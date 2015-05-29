//Ask user to choose the input and output directories
directory = getDirectory("Choose input directory");
fileList = getFileList(directory);
outputDirectory = getDirectory("Choose output directory");

//Ask user how many bins each stack should be divided into for each analysis
nBins = getNumber("The stacks should be divided into how many bins?", 5); 

//Set the image measurement tool to measure mean intensity for correct quantitation
run("Set Measurements...", "mean redirect=None decimal=9");

//Initialize result table
setResult("SampleID", 0, 0);
setResult("Number_of_Bins", 0, 0);
setResult("X_min_um", 0, 0);
setResult("X_max_um", 0, 0);
setResult("Y_min_um", 0, 0);
setResult("Y_max_um", 0, 0);
setResult("Z_min_um", 0, 0);
setResult("Z_max_um", 0, 0);

//Initialize arrays for storing mean intensity per bin results
XintensityArray = newArray(nBins);
YintensityArray = newArray(nBins);
ZintensityArray = newArray(nBins);

setBatchMode(true);

for (a=0; a<fileList.length; a++) {
	file = directory + fileList[a];

	//Open files in directory one at a time to be analyzed
	open(file);

	//++++++++++++Check data quality before processing+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	getVoxelSize(X, Y, Z, unit);
		
	//Double check to confirm voxels are cubic
	if(X/Y < 0.9 || Y/Z < 0.9 || Z/X < 0.9) {
		exit("Error: The voxel dimensions are not approximately cubic!");
	}
	
	//In first iteration, ask user to confirm voxel dimensions
	if(a == 0) {
		AverageDim = (X + Y + Z)/3;
		VoxelDim = getNumber("What are the cubic voxel dimensions in microns?", AverageDim);
	} 

	//Confirm that there are sufficient slices/pixels to have at least one slice/pixel per bin 
	getDimensions(width, height, channels, slices, frames);
	if (width < nBins || height < nBins || slices < nBins) {
			exit("Error: The # of bins is too large for this dataset!");
	}
	
	//++++++++++Binarize the CNT stack data for quantification+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//Set the max to a value of 1 and apply to whole stack
	setMinAndMax(0, 1); 
	run("Apply LUT", "stack");

	//Since apply rescales all values to 255, divide by 255 so that CNT voxel = 1, Polymer voxel = 0
	run("Divide...", "value=255 stack");

	//++++++++++++Measure percent intensity in X bins++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Generate a sum projection to allow for lateral quantitation of whole stack in one frame
	run("Z Project...", "projection=[Sum Slices]");

	//Divide slices into equal integer bins (round down to maintain equal # of slices in bins)
	XbinSize = floor(width/nBins);

	//Center the bins to be analyzed to middle of stack with unanalyzed remainder on top and bottom
	XstartPixel = floor((width%nBins)/2);
	
	//Calculate the Z positions of the start and end of the measure substack

	Xstart = XstartPixel * VoxelDim;
	Xend = (XstartPixel + nBins * XbinSize) * VoxelDim;
	
	for (b=0; b<nBins; b++) {
		//Calulate start point for iterative selection rectangle and create selection
		XbinStart = XstartPixel + XbinSize * b;		
		makeRectangle(XbinStart, 0, XbinSize, height);
	
		//measure mean in selection and insert mean result into Z array
		getStatistics(count, mean, min, max, std);
		XintensityArray[b] = mean/slices;
	}

	//++++++++++++Measure percent intensity in Y bins++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Divide slices into equal integer bins (round down to maintain equal # of slices in bins)
	YbinSize = floor(height/nBins);

	//Center the bins to be analyzed to middle of stack with unanalyzed remainder on top and bottom
	YstartPixel = floor((height%nBins)/2);
	
	//Calculate the Z positions of the start and end of the measure substack

	Ystart = YstartPixel * VoxelDim;
	Yend = (YstartPixel + nBins * YbinSize) * VoxelDim;
	
	for (b=0; b<nBins; b++) {
		//Calulate start point for iterative selection rectangle and create selection
		YbinStart = YstartPixel + YbinSize * b;		
		makeRectangle(0, YbinStart, width, YbinSize);
	
		//measure mean in selection and insert mean result into Z array
		getStatistics(count, mean, min, max, std);
		YintensityArray[b] = mean/slices;
	}

	//Close sum stack now that macro is finished analyzing it
	close("SUM*");


	//++++++++++++Measure percent intensity in Z bins++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//Divide slices into equal integer bins (round down to maintain equal # of slices in bins)
	ZbinSize = floor(nSlices/nBins);

	//Center the bins to be analyzed to middle of stack with unanalyzed remainder on top and bottom
	//Add 1 so for stack that is integer multiple of bins, start slice is 1 not 0
	startSlice = floor((nSlices%nBins)/2)+1;
	
	//Calculate the Z positions of the start and end of the measure substack

	Zstart = startSlice * VoxelDim;
	Zend = (startSlice + nBins * ZbinSize) * VoxelDim;
	
	for (b=0; b<nBins; b++) {

		//Reset the mean variable for each analysis iteration
		totalMean = 0;
		
		for (c=startSlice; c<(startSlice + ZbinSize); c++) {
			setSlice(c);
			getStatistics(count, mean, min, max, std);
			totalMean = mean + totalMean;
		}
		
		//Reset next starting slice to be one next slice in stack
		startSlice = c;
		
		//Calcualte mean Z intensity in bin and insert result into Z array
		ZintensityArray[b] = totalMean/ZbinSize;

	}

	//Close stack now that macro is finished analyzing it
	close();

	//++++++++++++++++++++++Output samples results to spreadsheet++++++++++++++++++++++++++++++++++++++++++++++
	//DO NOT DO THIS - this complicates finding the corresponding bundle csv file.  Only remove the extension, the R script can then directly edit
	//Get sample ID, removing extra contents from file name !!!Specific to this dataset!!!
	//SampleID = substring(fileList[a], 10, lengthOf(fileList[a])-12);
	
	setResult("SampleID", a, fileList[a]);	
	setResult("Number_of_Bins", a, nBins);
	setResult("X_min_um", a, Xstart);
	setResult("X_max_um", a, Xend);
	setResult("Y_min_um", a, Ystart);
	setResult("Y_max_um", a, Yend);
	setResult("Z_min_um", a, Zstart);
	setResult("Z_max_um", a, Zend);

	for (d=0; d<nBins; d++) {
		setResult("X_bin"+d+1, a, XintensityArray[d]);
	}
	
	for (d=0; d<nBins; d++) {
		setResult("Y_bin"+d+1, a, YintensityArray[d]);
	}
	
	for (d=0; d<nBins; d++) {
		setResult("Z_bin"+d+1, a, ZintensityArray[d]);
	}
}
setBatchMode(false);

//Save results table as csv file
if (nResults==0) exit("Results table is empty");
saveAs("Results",  outputDirectory + "CNT 3D Quantitation.csv"); 
run("Close");

