Note: 

To generate FigS2:

(1) Run [GESV_LRR.sh](https://github.com/peiyaohu/Hu_et_al_GESV_2025/blob/main/2_GESV_on_LRR-RLK/GESV_LRR.sh) to generate unique sequences of each sample

(2) The output of `GESV_LRR.sh` serves as the input file of [bwa.sh](https://github.com/peiyaohu/Hu_et_al_GESV_2025/blob/main/Figure1D-I/FigS2/bwa.sh). `bwa.sh` generates `.sam` files.

(3) The output of `bwa.sh` (`.sam` files) serves as the input file of [bwa_visual.R](https://github.com/peiyaohu/Hu_et_al_GESV_2025/blob/main/Figure1D-I/FigS2/bwa_visual.R).
