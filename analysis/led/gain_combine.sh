#! /bin/bash

# #data_file=20230428_0829_Calib
# #data_file=20230503_0737_Calib
# #data_file=20230503_0835_Calib
# #data_file=20230504_1425_Calib

# #data_file=20230516_1400_Calib
# #data_file=20230517_1550_Calib
# #data_file=20230517_1606_Calib
# #data_file=20230518_1533_Calib
# #data_file=20230519_1140_Calib
# #data_file=20230520_1556_Calib

# #beam=SPS
# beam=PS

# rm ../results/led/gain_allfilename.txt
# touch ../results/led/gain_allfilename.txt

# for txt in $( ls ../results/led/*_Calib_LED_all.root ); do
# #for txt in $( ls ../results/led/*5[12]*_Calib_LED_all.root ); do
# #  echo ${txt}
#   echo -n ${txt} >> ../results/led/gain_allfilename.txt
#   echo -n " " >> ../results/led/gain_allfilename.txt
# done

# var=`cat ../results/led/gain_allfilename.txt`
# ./between_files_gain ../results/led/gain_all.root ${var}



rm ../results/led/gain_allfilename.txt
touch ../results/led/gain_allfilename.txt

for txt in $( ls ../results/led/scan/*_Calib_LED.root ); do
#for txt in $( ls ../results/led/*5[12]*_Calib_LED_all.root ); do
#  echo ${txt}
  echo -n ${txt} >> ../results/led/gain_allfilename.txt
  echo -n " " >> ../results/led/gain_allfilename.txt
done

var=`cat ../results/led/gain_allfilename.txt`
./between_files_gain ../results/led/scan/all_days.root ${var}
