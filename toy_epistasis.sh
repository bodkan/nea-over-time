EPI_PAIRS=10
INIT_F=0.5
REPS=5

OUTPUT_DIR=simulations/toy_epistasis
mkdir -p $OUTPUT_DIR

for EPI_S in 1.0 0.1 0.05 0.01 0.005 0.001 0.0005 0.0001 0.0; do
    for RECOMB in 1e-6 1e-5 1e-4 1e-3; do
        for rep in `seq 1 $REPS`; do
            run_id="pairs_${EPI_PAIRS}__s_${EPI_S}__recomb_${RECOMB}__rep_${rep}"
            # qsub -V -cwd -j y -S /bin/bash \
     	    #      -l virtual_free=128MG,h_vmem=128M \
     	    #      -o ${tmp_dir}/sge/sge__${run_id}.out \
     	    #      -N ${run_id} \
 echo                 slim -d EPI_PAIRS=$EPI_PAIRS -d EPI_S=${EPI_S} -d RECOMB=$RECOMB -d INIT_F=$INIT_F \
                      -d "OUTPUT='${OUTPUT_DIR}/${run_id}.txt'" \
                      slim/toy_epistasis.slim
            exit
        done
    done
done

