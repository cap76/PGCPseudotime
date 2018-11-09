cat ../earls.bed ../lates.bed > PGC.bed
bedtools sort -i PGC.bed > PGCs.bed

computeMatrix scale-regions -S ../../MAX/SRR643235_trimmed_sorted.bw \
../../GSE79552/GSM2098095_NOMe-seq_Human_Week11_male_PGC_rep1.GCH.bw \
../../GSE79552/GSM2098097_NOMe-seq_Human_Week12_male_PGC_rep1.GCH.bw \
../../GSE79552/GSM2098099_NOMe-seq_Human_Week17_male_PGC_rep1.GCH.bw \
../../GSE79552/GSM2098105_NOMe-seq_Human_Week26_male_PGC_rep1.GCH.bw \
../../GSE79552/GSM2098095_NOMe-seq_Human_Week11_male_PGC_rep1.WCG.bw \
../../GSE79552/GSM2098097_NOMe-seq_Human_Week12_male_PGC_rep1.WCG.bw \
../../GSE79552/GSM2098099_NOMe-seq_Human_Week17_male_PGC_rep1.WCG.bw \
../../GSE79552/GSM2098105_NOMe-seq_Human_Week26_male_PGC_rep1.WCG.bw \
-R PGCs.bed --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --binSize 10 --outFileName MAXg --outFileNameMatrix MAXg_ --sortRegions no

computeMatrix scale-regions -S ../../MAX/SRR643235_trimmed_sorted.bw \
../../GSE79552/GSM2098095_NOMe-seq_Human_Week11_male_PGC_rep1.GCH.bw \
../../GSE79552/GSM2098097_NOMe-seq_Human_Week12_male_PGC_rep1.GCH.bw \
../../GSE79552/GSM2098099_NOMe-seq_Human_Week17_male_PGC_rep1.GCH.bw \
../../GSE79552/GSM2098105_NOMe-seq_Human_Week26_male_PGC_rep1.GCH.bw \
../../GSE79552/GSM2098095_NOMe-seq_Human_Week11_male_PGC_rep1.WCG.bw \
../../GSE79552/GSM2098097_NOMe-seq_Human_Week12_male_PGC_rep1.WCG.bw \
../../GSE79552/GSM2098099_NOMe-seq_Human_Week17_male_PGC_rep1.WCG.bw \
../../GSE79552/GSM2098105_NOMe-seq_Human_Week26_male_PGC_rep1.WCG.bw \
-R enhancers.bed --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --binSize 10 --outFileName MAXe --outFileNameMatrix MAXe_ --sortRegions no

computeMatrix scale-regions -S ../../MAX/SRR643235_trimmed_sorted.bw \
../../GSE79552/GSM2098095_NOMe-seq_Human_Week11_male_PGC_rep1.GCH.bw \
../../GSE79552/GSM2098097_NOMe-seq_Human_Week12_male_PGC_rep1.GCH.bw \
../../GSE79552/GSM2098099_NOMe-seq_Human_Week17_male_PGC_rep1.GCH.bw \
../../GSE79552/GSM2098105_NOMe-seq_Human_Week26_male_PGC_rep1.GCH.bw \
../../GSE79552/GSM2098095_NOMe-seq_Human_Week11_male_PGC_rep1.WCG.bw \
../../GSE79552/GSM2098097_NOMe-seq_Human_Week12_male_PGC_rep1.WCG.bw \
../../GSE79552/GSM2098099_NOMe-seq_Human_Week17_male_PGC_rep1.WCG.bw \
../../GSE79552/GSM2098105_NOMe-seq_Human_Week26_male_PGC_rep1.WCG.bw \
-R NPGCs.bed --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --binSize 10 --outFileName MAXnpcg --outFileNameMatrix MAXnpgc_ --sortRegions no

plotHeatmap -m MAXnpgc \
-out MAXnpgc.pdf \
--outFileSortedRegions MAXnpgc.bed  \
--interpolationMethod bilinear \
--whatToShow 'heatmap and colorbar' \
--samplesLabel MAX Wk11_GCH Wk12_GCH Wk17_GCH Wk26_GCH Wk11_WCG Wk12_WCG Wk17_WCG Wk26_WCG \
--colorMap Oranges \
--zMin 0 0 0 0 0 0 0 0 0 --zMax 13 .3 .3 .3 .3 .3 .3 .3 .3 \
--kmeans 4

plotProfile -m MAXnpgc -out MAXnpgc_.pdf --numPlotsPerRow 5  --samplesLabel MAX Wk11_GCH Wk12_GCH Wk17_GCH Wk26_GCH Wk11_WCG Wk12_WCG Wk17_WCG Wk26_WCG \
--yMin 0 0 0 0 0 0 0 0 0 --yMax 13 .3 .3 .3 .3 .3 .3 .3 .3 \

#plotHeatmap -m MAXg \
#-out MAXg.pdf \
#--outFileSortedRegions MAXg.bed  \
#--interpolationMethod bilinear \
#--whatToShow 'heatmap and colorbar' \
#--samplesLabel MAX Wk11_GCH Wk12_GCH Wk17_GCH Wk26_GCH Wk11_WCG Wk12_WCG Wk17_WCG Wk26_WCG \
#--colorMap Oranges \
#--zMin 0 0 0 0 0 0 0 0 0 --zMax 13 .3 .3 .3 .3 .3 .3 .3 .3 \
#--kmeans 4

#plotProfile -m MAXg -out MAXg_.pdf --numPlotsPerRow 5  --samplesLabel MAX Wk11_GCH Wk12_GCH Wk17_GCH Wk26_GCH Wk11_WCG Wk12_WCG Wk17_WCG Wk26_WCG \
#--yMin 0 0 0 0 0 0 0 0 0 --yMax 13 .3 .3 .3 .3 .3 .3 .3 .3 \

