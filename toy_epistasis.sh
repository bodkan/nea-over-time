EPI_PAIRS=1000
INIT_F=0.03
REPS=20

OUTPUT_DIR=simulations/toy_epistasis
mkdir -p $OUTPUT_DIR

for EPI_S in 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2; do
    for EPI_REC in 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2; do
        for rep in `seq 1 $REPS`; do
            run_id="n_${EPI_PAIRS}__s_${EPI_S}__r_${EPI_REC}__f_${INIT_F}__rep_${rep}"
            qsub -V -cwd -j y -S /bin/bash \
     	         -l virtual_free=256M,h_vmem=256M \
     	         -o tmp/sge/sge__${run_id}.out \
     	         -N ${run_id} \
                 -b y \
                 `which slim` -d EPI_PAIRS=$EPI_PAIRS -d EPI_S=${EPI_S} -d EPI_REC=$EPI_REC -d INIT_F=$INIT_F \
                      -d "OUTPUT='${OUTPUT_DIR}/${run_id}.txt'" \
                      slim/toy_epistasis.slim
        done
    done
done
