cd D:/Research/ogrdbstats/testdata/VH_tigger
Rscript D:/Research/ogrdbstats/testdata/ogrdbstats.R IMGT_REF_GAPPED.fasta Homosapiens TWO01A_naive_genotyped.tsv VH --inf TWO01A_naive_novel_ungapped.fasta --hap IGHJ6

cd D:/Research/ogrdbstats/testdata/V_tigger_truncated
Rscript D:/Research/ogrdbstats/testdata/ogrdbstats.R IMGT_REF_GAPPED.fasta Homosapiens  TWO01A_naive_genotyped.tsv VH  --inf TWO01A_naive_novel_ungapped.fasta --hap IGHJ6

cd D:/Research/ogrdbstats/testdata/JH_tigger
Rscript D:/Research/ogrdbstats/testdata/ogrdbstats.R IMGT_REF_GAPPED_J_CHANGES.fasta Homosapiens TWO01A_naive.airr.tab JH  --inf TWO01A_naive_novel.fasta --hap IGHV2-5

cd D:/Research/ogrdbstats/testdata/D_tigger
Rscript D:/Research/ogrdbstats/testdata/ogrdbstats.R IMGT_REF_GAPPED_D_1_26_01_removed.fasta Homosapiens TWO01A_naive.airr.tab DH --inf TWO01A_naive_novel.fasta --hap IGHJ6

cd D:/Research/ogrdbstats/testdata/JH_igdiscover
Rscript D:/Research/ogrdbstats/testdata/ogrdbstats.R IMGT_REF_GAPPED_fake_j.fasta Homosapiens filtered.tab JH --inf J.fasta  --hap IGHV2-5

cd D:/Research/ogrdbstats/testdata/JL_igdiscover
Rscript D:/Research/ogrdbstats/testdata/ogrdbstats.R IMGT_REF_GAPPED_fake_j_fake_JK.fasta Homosapiens filtered.tab JH --inf J.fasta  --hap IGHV2-5

cd D:/Research/ogrdbstats/testdata/VH_partis
Rscript D:/Research/ogrdbstats/testdata/ogrdbstats.R IMGT_REF_GAPPED.fasta Homosapiens TW02A_OGRDB.tsv VH --inf TW02A_V_OGRDB.fasta --hap IGHJ6
