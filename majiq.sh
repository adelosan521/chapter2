##This code uses Majiq to build and analyze GFF3 files from Wilfried Haerty and then uses Voila to view and dump the results to a tsv file.

##Majiq
##Used GFF3 from Wilfried Haerty ("untar" the following file: Talon_minID90_N5K3_purged_filtered_gencode_renamed_validated_NoSequins.gff.tar.gz)

##Build the MAJIQ splice graph
majiq build /t1-data/user/aangeles/fastaq_files/MAJIQ/Talon_minID90_N5K3_purged_filtered_gencode_renamed_validated_NoSequins.gff -o /t1-data/user/aangeles/fastaq_files/MAJIQ/build -c /t1-data/user/aangeles/fastaq_files/MAJIQ/majiq_build_config4.ini

##Calculate the delta PSI (percent spliced in) values between the two gruops of samples
majiq deltapsi -o /t1-data/user/aangeles/fastaq_files/MAJIQ/deltapsi -grp1 /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914279_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914288_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914293_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914305_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914306_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914307_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914316_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914317_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914321_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914322_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914325_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914326_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914328_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914330_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914332_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914334_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914335_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914336_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914345_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914346_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR914347_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR946960_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR946971_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR946973_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR946974_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR946976_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR946978_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR946980_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR946987_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR947002_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR947004_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR947006_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR947007_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR947008_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR947013_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR947014_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR947016_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR947017_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1243458_Aligned.out.majiq -grp2  /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775544_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775551_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775552_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775553_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775554_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775555_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775556_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775562_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775563_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775564_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775565_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775566_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775591_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775592_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775593_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775595_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775596_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775597_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775598_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775599_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775600_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775601_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775602_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775631_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775634_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775636_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775637_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775638_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775639_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775640_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775641_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775643_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775644_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775685_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775686_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775687_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775688_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775689_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775692_Aligned.out.majiq /t1-data/user/aangeles/fastaq_files/MAJIQ/build/ERR1775693_Aligned.out.majiq  -n ipsc neuron

##Visualize the results using Voila (creates visual representation of splicing events using Voila viewer)
voila view /t1-data/user/aangeles/fastaq_files/MAJIQ/deltapsi

##Export Voila results to a TSV file 
voila tsv -f /t1-data/user/aangeles/fastaq_files/MAJIQ/build/voila_output.tsv /t1-data/user/aangeles/fastaq_files/MAJIQ/build/


##majiq deltapsi command from ceph
majiq deltapsi --min-deltapsi 0.05 -o /t1-data/user/aangeles/fastaq_files/MAJIQ/deltapsi -grp1 /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914279_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914288_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914293_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914305_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914306_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914307_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914316_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914317_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914321_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914322_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914325_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914326_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914328_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914330_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914332_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914334_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914335_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914336_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914345_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914346_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR914347_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR946960_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR946971_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR946973_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR946974_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR946976_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR946978_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR946980_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR946987_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR947002_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR947004_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR947006_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR947007_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR947008_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR947013_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR947014_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR947016_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR947017_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR1243458_Aligned.out.majiq -grp2  /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR1775544_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR1775551_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR1775552_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR1775553_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR1775554_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR1775555_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR1775556_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR1775562_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR1775563_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775564_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775565_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775566_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775591_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775592_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775593_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775595_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775596_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775597_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775598_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775599_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775600_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775601_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775602_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775631_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775634_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775636_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775637_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775638_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775639_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775640_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775641_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775643_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775644_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775685_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775686_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775687_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR1775688_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775689_Aligned.out.majiq /ceph/project/tunbridgelabaangeles/fastaq_files/MAJIQ/build/ERR1775692_Aligned.out.majiq /ceph/project/tunbridgelab/aangeles/fastaq_files/MAJIQ/build/ERR1775693_Aligned.out.majiq  -n ipsc neuron


xx
