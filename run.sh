echo "Running data generation ..." 


date_time="report_"
date_time=$date_time`date +"%Y_%m_%d_%H_%M_%S"`


cd r_code
echo "Creating report folder: $date_time ..."
mkdir $date_time


echo "Running classification ..."
echo "R script output:"
echo "-----------------------------------------------------------------------"
Rscript "main.R" "donor.dat" $date_time "sv"
echo "-----------------------------------------------------------------------"


cd ..
echo "Done!"


