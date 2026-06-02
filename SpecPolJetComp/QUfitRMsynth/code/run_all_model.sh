# Note WAPITI/QU names just map onto different computers used to process the data

# Bright models for the bulk of the flaring epochs
for MODEL in S SS SSS T ST P PP;do
    echo "Running model $MODEL..."
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_WAPITI_20230904.txt
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_QU_20230906.txt
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_QU_20230908.txt
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_QU_20230916.txt
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_WAPITI_20230923.txt
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_QU_20231001.txt
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_WAPITI_20231016.txt
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_QU_20231022.txt
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_QU_20231028.txt
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_QU_20231106.txt
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_QU_20231112.txt
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_QU_20231118.txt
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_QU_20231125.txt
done

# Peak 1 (fainter peak) has extra models
for MODEL in S SS SSS SSSS SSSSS T ST SST STT SSST SSTT STTT P PP PPP;do
    echo "Running model $MODEL for peak 1..."
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_WAPITI_20231006.txt
done

# Peak 2 (brighter peak) has extra models
for MODEL in S SS SSS SSSS SSSSS SSSSSS T ST SST STT SSST SSTT STTT SSSST SSSTT SSTTT STTTT P PP PPP PPPP;do
    echo "Running model $MODEL for peak 2..."
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_WAPITI_20231014.txt
done

# Post-flaring epochs (core)
for MODEL in S SS T ST P;do
    echo "Running model $MODEL for post-flaring epochs..."
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_QU_20240210.txt
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_QU_20240219.txt
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_QU_20240225.txt
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_QU_20240331.txt
done

# Post-flaring epochs (ejecta)
for MODEL in S T P;do
    echo "Running model $MODEL for post-flaring epochs (ejecta)..."
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_QU_20240210_ejecta.txt
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_QU_20240219_ejecta.txt
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_QU_20240225_ejecta.txt 
    python3 faraday_obs_fit.py $MODEL ../files/J1727/QU_text/SwiftJ1727_QU_20240331_ejecta.txt
done