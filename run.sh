# Skrypt uruchomieniowy opakowujący wywołania.
# Autorzy: Patryk Pankiewicz, Łukasz Brzezicki


# get report dir name
report_dir="report_"
report_dir=$report_dir`date +"%Y_%m_%d_%H_%M_%S"`

# test data
data="data/araclean.dat"

#result file 
seq="donor.dat"

# A length
A=50

# B length
B=5

# type
type="ACCEPTOR"

# additional report file
report_file="report.txt"

# overlap
overlap="-o"

# validation method
validation="cv"

echo "Running data generation ..." 

python "python_code/get_acceptors_and_donors.py" $overlap -A $A -B $B -i $data -t $type -r $seq
 
echo "Creating report folder: $date_time ..."
mkdir $report_dir

cd "./r_code"

echo "Running classification ..."
echo "R script output:"
echo "-----------------------------------------------------------------------"
Rscript "main.R" "../"$seq "../"$report_dir $validation
echo "-----------------------------------------------------------------------"

cd ..

echo "A length: $A, B length: $B" >> $report_dir"/"$report_file   
echo "Type: $type" >> $report_dir"/"$report_file   
echo "Overlap: $overlap" >> $report_dir"/"$report_file
echo "Validation: $validation" >> $report_dir"/"$report_file

mv $seq $report_dir"/"

echo "Done!"


