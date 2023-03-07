cd VH_airr_no_novel
Rscript ../../ogrdbstats.R IMGT_REF_GAPPED.fasta Homosapiens P8_I1_S1_airr.tsv VH
cd..

cd VH_tigger
Rscript ../../ogrdbstats.R IMGT_REF_GAPPED.fasta Homosapiens TWO01A_naive_genotyped.tsv VH --inf_file TWO01A_naive_novel_ungapped.fasta --hap_gene IGHJ6
cd ..

cd V_tigger_truncated
Rscript ../../ogrdbstats.R IMGT_REF_GAPPED.fasta Homosapiens  TWO01A_naive_genotyped.tsv VH  --inf_file TWO01A_naive_novel_ungapped.fasta --hap_gene IGHJ6
cd ..

cd JH_tigger
Rscript ../../ogrdbstats.R IMGT_REF_GAPPED_J_CHANGES.fasta Homosapiens TWO01A_naive.airr.tab JH  --inf_file TWO01A_naive_novel.fasta --hap_gene IGHV2-5
cd ..

cd D_tigger
Rscript ../../ogrdbstats.R IMGT_REF_GAPPED_D_1_26_01_removed.fasta Homosapiens TWO01A_naive.airr.tab DH --inf_file TWO01A_naive_novel.fasta --hap_gene IGHJ6
cd ..

cd JH_igdiscover
Rscript ../../ogrdbstats.R IMGT_REF_GAPPED_fake_j.fasta Homosapiens filtered.tab JH --inf_file J.fasta  --hap_gene IGHV2-5
cd ..

cd JL_igdiscover
Rscript ../../ogrdbstats.R IMGT_REF_GAPPED_fake_j_fake_JK.fasta Homosapiens filtered.tab JK --inf_file J.fasta  --hap_gene IGHV2-5
cd ..

cd VH_partis
Rscript ../../ogrdbstats.R IMGT_REF_GAPPED.fasta Homosapiens TW02A_OGRDB.tsv VH --inf_file TW02A_V_OGRDB.fasta --hap_gene IGHJ6
cd ..

cd "private/PRJEB30386 - Kappa"
Rscript ../../../ogrdbstats.R IMGT_REF_GAPPED.fasta Homosapiens Read_file.tab VK --inf_file Inferred_file.fasta --hap_gene IGHJ6
cd ..
cd ..

cd "private/hamster" 
Rscript ../../../ogrdbstats.R hamster_IGH_VDJ.fasta Hamster igblast_001_x-clones.tsv VH --all_novel
cd ..
cd ..