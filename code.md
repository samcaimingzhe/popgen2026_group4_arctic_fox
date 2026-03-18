Copy the data to your path.
```
cp -r /course/popgenmsc26/projects/arctic_fox/data/* ./
```
Change the name of chromosomes into integer with python script.
```python
with open("AF.imputed.thin.0.bim", "r") as infile:
    content = infile.read()

for i in range(1, 25):
    acc = f"NC_{54823+i:06d}.1"
    content = content.replace(acc, str(i))

with open("AF.imputed.thin.bim", "w") as outfile:
    outfile.write(content)

print("Done!")
```
```
cp AF.imputed.thin.bim AF.imputed.thin.0.bim
python3 chr.py
```

Obtain PCA results by plink. 10 PCs will be calculted though only first 2 are useful.
```
plink --bfile AF.imputed.thin --freq --chr-set 24 --out AF_freqs
plink --bfile AF.imputed.thin --read-freq AF_freqs.frq --pca 10 --chr-set 24 --out AF_pca
```
Example code to download the pca data files to your local for plotting.

NB: remember to change the user number!
```
mkdir -p ~/Desktop/AF
scp "popgenmsc26user00@emily.popgen.dk:~/github/AF_pca.*" ~/Desktop/AF
scp "popgenmsc26user00@emily.popgen.dk:~/github/sample_popinfo.tsv" ~/Desktop/AF
```
```R
setwd('~/Desktop/AF') # make a AF folder on your desktop
data = read.table('AF_pca.eigenvec')
popinfo = read.table('sample_popinfo.tsv',header = T)
colnames(data) <- c('fid', 'iid', 
                    'pc1', 'pc2', 'pc3','pc4', 'pc5',
                    'pc6','pc7', 'pc8', 'pc9','pc10')
data$fid = popinfo$Region

library(plotly)
plot_ly(data, 
        x = ~pc1, 
        y = ~pc2, 
        z = ~pc3, 
        color = ~fid,     
        colors = "Set1",  
        type = 'scatter3d', 
        mode = 'markers',
        marker = list(size = 4, opacity = 0.8),
        text = ~iid)
```
<img width="729" height="476" alt="pca" src="https://github.com/user-attachments/assets/82f2f873-3922-41ac-bc99-f4803806e810" />
So we observed that ArcFoxSamp_128, ArcFoxSamp_144, ArcFoxSamp_152 may be the 3 most recent immigrants.

Admixture assume the SNPs are independent to each others, so we have to do LD prunning.
```
# remove SNPs with missing genotypes
plink --bfile AF.imputed.thin --geno 0 --make-bed --out AF.complete
# remove low MAF SNPs
plink --bfile AF.complete --maf 0.05 --make-bed --out AF.frq_normal
# LD pruning
plink --bfile AF.frq_normal --indep-pairwise 50 10 0.1 --out AF.LD
plink --bfile AF.frq_normal --extract AF.LD.prune.in  --make-bed --out AF.clean
```
Then we can check how many SNPs are filtered out during each step.
```
wc -l *.bim|sort -r
```
```
 2906282 total
  816284 AF.imputed.thin.bim
  816284 AF.imputed.thin.0.bim
  719391 AF.complete.bim
  498777 AF.frq_normal.bim
   55546 AF.clean.bim
```
Then we can run the admixture by K loop and run loop, try to run it in `screen`:
```
screen -S admix
```
Inside admix, we run this double loop commands:
```
for K in {2..7}; do
  for i in {1..100}; do
    echo "running K=${K} & run=${i}"
    admixture -s ${i} AF.clean.bed ${K} > AF.clean_K${K}_run${i}.log
    cp AF.clean.${K}.Q AF.clean_K${K}_run${i}.Q
    cp AF.clean.${K}.P AF.clean_K${K}_run${i}.P          
  done
done
```
You may use `control + A + D` to exit the screen and log out the server, because admixture requires a while to run.
We can continue to run Fst test. Please download the files to local.
```
scp "popgenmsc26user00@emily.popgen.dk:~/github/AF.imputed.thin.???" ~/Desktop/AF
```






