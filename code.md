Copy the data to your path
cp -r /course/popgenmsc26/projects/arctic_fox/data/* ./

Obtain PCA results by plink

plink --bfile AF.imputed.thin --freq --chr-set 24 --out my_freqs
plink --bfile AF.imputed.thin --read-freq my_freqs.frq --pca 3 --chr-set 24








