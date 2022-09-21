#!/bin/bash
result=$(wc -l < run_batch_simulations.sh)
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 29 0.0 0.1 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 29 0.0 1 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 29 0.0 10 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 29 0.00001 0.0001 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 29 0.00001 0.0001 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 29 1e-5 0.001 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 29 1e-5 0.01 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 29 1e-5 0.0 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 29 1e-5 0.1 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 29 1e-5 1 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 29 1e-5 10 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 1 0 29 29 0 0 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 1 0 29 29 1e-5 0 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 1 0 29 29 1e-5 0.1 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 1 0 29 29 1e-5 1 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 1 0 29 29 1e-5 10 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 12 0 29 29 0 0 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 12 0 29 29 1e-5 0 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 12 0 29 29 1e-5 0.1 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 12 0 29 29 1e-5 1 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 12 0 29 29 1e-5 10 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 6 0 29 29 0 0 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 6 0 29 29 1e-5 0 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 6 0 29 29 1e-5 0.1 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 6 0 29 29 1e-5 1 1 306 12
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 6 0 29 29 1e-5 10 1 306 12
echo "Line: ${LINENO} out of ${result}"

bash make_df_simulations.sh
bash make_batch_movies.sh

#python test model_active.py k_g1 k_s k_g2 k_m D a run duration sim_type