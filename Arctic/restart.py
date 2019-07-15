import subprocess
import os.path
import netCDF4

num = input("Restart number? ")
file_path = "$AKDRIVE/Arctic4/run16/restart_" + num
if os.path.exists(file_path):
    print("directory exists already:", file_path)
    exit()
else:
    cmd = "mkdir " + file_path
    subprocess.call([cmd], shell=True)

cmd = "mv -i arctic4_avg_* tmp"
subprocess.call([cmd], shell=True)
cmd = "mv -i arctic4_sta.nc arctic4_sta_" + num + ".nc"
subprocess.call([cmd], shell=True)
#fh = netCDF4.Dataset("arctic4_flt.nc", "r")
#nt = len(fh.dimensions['ocean_time'])
#nt = str(nt-1)
#print(nt)
#fh.close()
#cmd = "mv -i arctic4_flt.nc arctic4_flt_" + num + ".nc"
#subprocess.call([cmd], shell=True)
#cmd = "ncks -d ocean_time,"+nt+","+nt+" arctic4_flt_" + num + ".nc arctic4_flt.nc"
#subprocess.call([cmd], shell=True)
cmd = "cp *rst* " + file_path
subprocess.call([cmd], shell=True)
#vi ocean_arctic4.in
#cp *avg_* $AKDRIVE/Arctic4/run16/averages
#cp *avg2_* $AKDRIVE/Arctic4/run16/averages2
