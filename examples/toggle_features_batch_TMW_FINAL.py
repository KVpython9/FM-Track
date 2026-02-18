"""
    toggle_features.py: an example demonstrating how to run FM-Track. It runs on the example
    data placed inside of the examples folder of FM-Track. This file demonstrates
    how to modify many important parameters in FM-Track.
"""


import fmtrack
import numpy as np
import os
import datetime
import traceback
import glob


path =[r"/Volumes/SG - 500G/mMVICS CytoD/20260209n1 mMVICS WT P5 CytoD 1PT03_BS-PT863_VS_Bconc_0PT1/Cell 22/",
       #r"/Volumes/SG - 500G/mMVICS CytoD/20260209n1 mMVICS WT P5 CytoD 1PT03_BS-PT863_VS_Bconc_0PT1/Cell 2/",
       #r"/Volumes/SG - 500G/mMVICS CytoD/20260209n1 mMVICS WT P5 CytoD 1PT03_BS-PT863_VS_Bconc_0PT1/Cell 3/",
       #r"/Volumes/SG - 500G/mMVICS CytoD/20260209n1 mMVICS WT P5 CytoD 1PT03_BS-PT863_VS_Bconc_0PT1/Cell 4/",
       #r"/Volumes/SG - 500G/mMVICS CytoD/20260209n1 mMVICS WT P5 CytoD 1PT03_BS-PT863_VS_Bconc_0PT1/Cell 5/",
       #r"/Volumes/SG - 500G/mMVICS CytoD/20251211n1 mMVICS WT P5 CytoD 1PT03_BS-PT863_VS/Cell 2/",
       #r"/Volumes/SG - 500G/mMVICS CytoD/20251211n1 mMVICS WT P5 CytoD 1PT03_BS-PT863_VS/Cell 3/",
       #r"/Volumes/SG - 500G/mMVICS CytoD/20251211n1 mMVICS WT P5 CytoD 1PT03_BS-PT863_VS/Cell 4/",
       #r"/Volumes/SG - 500G/mMVICS CytoD/20251211n1 mMVICS WT P5 CytoD 1PT03_BS-PT863_VS/Cell 5/",
       #r"/Volumes/SG - 500G/mMVICS CytoD/20251211n1 mMVICS WT P5 CytoD 1PT03_BS-PT863_VS/Cell 6/"
       #r"/Volumes/SG - 500G/mMVICS CytoD/20251211n1 mMVICS WT P5 CytoD 1PT03_BS-PT863_VS/Cell 7/",
       #r"/Volumes/SG - 500G/mMVICS CytoD/20251211n1 mMVICS WT P5 CytoD 1PT03_BS-PT863_VS/Cell 8/",
       #r"/Volumes/SG - 500G/mMVICS CytoD/20251211n1 mMVICS WT P5 CytoD 1PT03_BS-PT863_VS/Cell 9/",
       #r"/Volumes/SG - 500G/mMVICS CytoD/20251211n1 mMVICS WT P5 CytoD 1PT03_BS-PT863_VS/Cell 10/",
       #r"/Volumes/SG - 500G/mMVICS CytoD/20251211n1 mMVICS WT P5 CytoD 1PT03_BS-PT863_VS/Cell 11/",
       #r"/Volumes/SG - 500G/mMVICS CytoD/20251211n1 mMVICS WT P5 CytoD 1PT03_BS-PT863_VS/Cell 12/",
       ]

should_plot = False   # toggles on and off 3D plotting
        
# set microscope field-of-view dimensions
# core Zeiss LSM 710 is 149.95, 149.95, 140
# Zoldan Leica LSM is 174.94, 174.94, 91.2 or 100.8
# core Nikon NSPARC single cell is 147, 147, 119.23 (note that the setting is 220.58 but the spherical aberration is 1.85 so 220.58/1.85 = 119.23; this gives us isotropic voxels)
# core Nikon NSPARC networks is 441.856, 441.856, 119.23 (same note as above)

dims = np.array([441.856, 441.856, 119.23])

          

for j in range(len(path)):
    try:
        # (2) find bead positions from images
        print('Importing bead data for j{} at {}.' .format(j, datetime.datetime.now().strftime("%H:%M")))
        bead_channel = 1
        beads_init = fmtrack.FMBeads()
        beads_init.get_bead_centers(path[j] + "Drug/Beads/Drug.nd2 - Drug%s.tif", dims)
        #bead_radius = 0.5          # in pixels (CHANGE THIS) added by kv
        #intensity_threshold = 300 #added by kv
        
        
        beads_final = fmtrack.FMBeads()
        beads_final.get_bead_centers(path[j] + "Basal/Beads/Basal.nd2 - Basal%s.tif", dims)
        #bead_radius=0.5     #added by kv
        #intensity_threshold=300    #added by kv
        


        # (4) run tracking algorithm
        tracker = fmtrack.FMTracker(beads_init=beads_init, beads_final=beads_final)
        tracker.print_progress = True   # toggles on and off printing status bars
        tracker.track_type = 1   # 1 = no tracking correction, 2 = tracking correctiion for microscope drift
        tracker.num_feat = 5   # number of feature vectors
        tracker.num_nearest = 5   # number of beads in final state to compare to a particular bead in the initial state 
        tracker.use_box = False   # True = far field boundary is rectangular prism, False = far field boundary is determined by dist. to cell surface
        tracker.run_tracking()
        
        # (5) save all of the output data
        track_folder = path[j] + "Drug/Cells/Smoothed_SA_Corrected/Track"
        tracker.save_all(track_folder, dims=dims)
        
        print("Track for j{} complete and saved at {}." .format(j, datetime.datetime.now().strftime("%H:%M")))
        
        ######################################################################################################
    except Exception as e:
        # Handle errors if needed
        date = datetime.datetime.now().strftime("%Y%m%d")
        print("Skipping j{} due to an error at {}. Continuing onto j{}." .format(j, date, j+1))

        error_text = traceback.format_exc()
        error_file_path = os.path.join(path[j], "error_log_{}_for_j{}.txt" .format(date, j))
        with open(error_file_path, "w") as error_file:
            error_file.write(error_text)

        print("Error details saved to {}, which reads as follows: {}" .format(j, e))
        print("")
        print("Program will continue with other files.")
        print("")


print("All files complete.")
