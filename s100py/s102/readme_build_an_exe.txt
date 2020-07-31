create a nomkl environment (mkl takes hundreds of megs compressed -- far more than the program warrants)
copy the s100py folder into the lib/site-packages
then use pyinstaller to make an executable

conda create -n nomkl python=3.6 nomkl pyinstaller gdal h5py numpy

Now copy s100py to envs/nomkl/lib/site-packages/s100py and 
change directroy to envs/nomkl/lib/site-packages/s100py

pyinstaller --onefile -n make_s102 utils.py 
(if proj.db isn't being found then add the data and modify the utils.py script to find it)
(to debug if the files are included use --onedir instead of --onefile)
pyinstaller --onefile --add-data C:\PydroXL_19_Dev\envs\nomkl2\Library\share\proj;Library\share\proj -n make_s102 utils.py

The resulting executable will be at dist/make_s102.exe