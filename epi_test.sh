EPI_PAIRS=1000
INIT_F=0.5

for EPI_S in 1.0 0.5 0.1 0.05 0.01 0.005 0.001 0.0005 0.0001 0.0; do
    for RECOMB in 1e-6 1e-5 1e-4 1e-3; do
        slim -d EPI_PAIRS=$EPI_PAIRS -d EPI_S=-$EPI_S -d RECOMB=$RECOMB -d INIT_F=$INIT_F slim/toy_epistasis.slim
        exit
    done
done

