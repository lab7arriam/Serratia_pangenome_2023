#!/usr/bin/zsh

mkdir -p serratia_refs_190921
ftps=( $(esearch -db assembly -query '"Serratia marcescens"[Organism] AND (latest[filter] AND "complete genome"[filter] AND all[filter] NOT anomalous[filter] AND "refseq has annotation"[Properties])' | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq) )
for ftp in ${ftps[@]}; do
    acc=$(echo ${ftp} | rev | cut -d'/' -f1 | rev )
    mkdir -p serratia_refs_190921/${acc}
    wget -P serratia_refs_190921/${acc} ${ftp}/${acc}_genomic.fna.gz
    wget -P serratia_refs_190921/${acc} ${ftp}/${acc}_genomic.gff.gz
done
