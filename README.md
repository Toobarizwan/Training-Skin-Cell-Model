# Training-Skin-Cell-Model
A skin-cell-specific TWAS pipeline where I train ctPred models on normal human skin scRNA-seq data (800k+ cells, 34 datasets)

Installing Pipeline:
````
git clone https://github.com/Toobarizwan/Training-Skin-Cell-Model
````


## Run NEXTFLOW
This step will produce .db files.

````
conda activate scPrediXcan
````

````
setsid nohup bash /home/tooba/Skin_Cell_Model/NEXTFLOW_scripts/Venous_EC_nextflow.sh > /home/tooba/Skin_Cell_Model/predictdb.logs/Venous_EC_nextflow.log 2>&1 < /dev/null &
disown
````
