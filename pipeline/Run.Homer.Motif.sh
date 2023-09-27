#!/bin/sh

findMotifsGenome=/mnt/d/WSL2/homer/bin/findMotifsGenome.pl
bed_dir=/mnt/e/LabWork/Project/CRC_NGS_ATAC/CRC_Epi_scATAC/Results/05.Epi_TF_Clustering/diff_peak_cluster/bed

ls ${bed_dir} |
    awk '{gsub(".bed", "");gsub("marker.peak.tumor.", ""); print $1}' > sample.list
split -l 4 sample.list -d -a 1  sample.list

for i in {0..6}
do  
    echo -e "["$(date)"]\tStart job submission round ${i}...." 

    while read sample; do
        nohup ${findMotifsGenome} \
            ${bed_dir}/marker.peak.tumor.${sample}.bed \
            hg38 homer/clusters/${sample} \
            -size 200 \
            > homer/homer.${sample}.log 2>&1 &

        echo -e "["$(date)"]\tSubmit job for ${sample}...."
done < sample.list${i};
sleep 60m
done

