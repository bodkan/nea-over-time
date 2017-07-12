EPI_PAIRS=1000
INIT_F=0.5
REPS=20

OUTPUT_DIR=simulations/toy_epistasis
mkdir -p $OUTPUT_DIR

for EPI_S in 0.01 0.005 0.001 0.0005 0.0001 0.00005 0.00001; do
    for EPI_REC in 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3; do
        for rep in `seq 1 $REPS`; do
            run_id="n_${EPI_PAIRS}__s_${EPI_S}__r_${EPI_REC}__f_${INIT_F}__rep_${rep}"
            qsub -V -cwd -j y -S /bin/bash \
     	         -l virtual_free=128M,h_vmem=128M \
     	         -o tmp/sge/sge__${run_id}.out \
     	         -N ${run_id} \
                 -b y \
                 `which slim` -d EPI_PAIRS=$EPI_PAIRS -d EPI_S=${EPI_S} -d EPI_REC=$EPI_REC -d INIT_F=$INIT_F \
                      -d "OUTPUT='${OUTPUT_DIR}/${run_id}.txt'" \
                      slim/toy_epistasis.slim
        done
    done
done
