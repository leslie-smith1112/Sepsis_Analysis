filenames=`ls ./*/metadata_*`
target="trimmed_metadata"
for eachfile in $filenames
do
cut -f1,2,6,7,10,11,16,18,19,27 $eachfile > "$eachfile$target"
done
