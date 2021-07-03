# Fiesta_motility_analysis
Processing and analysis of tracking data from FIESTA :

These codes use tracking data from FIESTA. Individual tracks should be saved as .mat files in the designated folder (see line 22 for change_point_analysis and line 38 for reversal_event_analysis).

Use still images of polarity marked microtubules to determine the direction of motility. Set the scale of the image to microns. For each track use the multi-point tool in ImageJ to mark the coordinates of the microtubule. The first point should be the minus-end of the microtubule. Save the coordinates as a .xlsx file in the same folder.

Run ‘change_point_analysis.m’ to segment tracks by change points. 
Lines 59-64 are input parameters that may need to be adjusted to run the change point function. 
This code runs tmsd_changepts function to find local alpha-values using a rolling MSD over a window. The tmsd_changepts function will output a matrix with x, y, alpha values, and categorize individual runs as processive (1), or diffusive (0). The change_point_analysis code will use this information to find processive and diffusive runs and perform further analysis. 

Run ‘reversal_event_analysis.m’ to segment tracks by reversals. 
Lines 277 and 278 are segments of the code that find stationary, diffusive, and processive periods of motility based on the absolute value of the run-length in between two reversals. Currently, the code finds:
Stationary runs = run length <=10 nm,
Diffusive runs = 200nm > run length > 10 nm,
Processive runs = run length > 200 nm.
These values can be adjusted for different experimental set-ups.

Included in the folder are two timelapse recordings and tracking data of cargoes moving along microtubules; ‘phago_minus_2’ is a cargo that moves uni-directionally towards the minus-end and ‘phago_bidir_1’ is a cargo that moves bi-directionally. The microtubule coordinates are saved as .xlsx files for each example. 
