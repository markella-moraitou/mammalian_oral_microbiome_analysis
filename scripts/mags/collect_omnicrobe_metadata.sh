##### Collect metadata from omnicrobe database #####

#### Use taxids to get info on known habitats, phenotypes and use for all assembled MAGs

#### VARIABLES AND WORKING DIRECTORY ####
scriptdir=".."
outdir="../../output/mags"

#### PREP INPUT ####
cut -d, -f 3 $outdir/name_to_taxid.csv > $outdir/taxids.tmp
# Remove NAs and header
sed -i '/NA/d' $outdir/taxids.tmp
sed -i '/taxids/d' $outdir/taxids.tmp

#### RUN SCRIPT ####
python $scriptdir/access_omnicrobe_db.py $outdir/taxids.tmp $outdir "habitat"
python $scriptdir/access_omnicrobe_db.py $outdir/taxids.tmp $outdir "phenotype"
python $scriptdir/access_omnicrobe_db.py $outdir/taxids.tmp $outdir "use"

rm $outdir/taxids.tmp
