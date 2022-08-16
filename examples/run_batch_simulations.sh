#!/bin/bash
result=$(wc -l < run_batch_simulations.sh)
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 1 0 29 29 0 0 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 6 0 29 29 0 0 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 12 0 29 29 0 0 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 29 0 0.1 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 29 0 1 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 29 0 10 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 1 0 29 29 0.00001 0 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 1 0 29 29 0.00001 0.1 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 1 0 29 29 0.00001 1 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 1 0 29 29 0.00001 10 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 6 0 29 29 0.00001 0 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 6 0 29 29 0.00001 0.1 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 6 0 29 29 0.00001 1 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 6 0 29 29 0.00001 10 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 12 0 29 29 0.00001 0 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 12 0 29 29 0.00001 0.1 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 12 0 29 29 0.00001 1 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 12 0 29 29 0.00001 10 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 29 0.00001 0 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 29 0.00001 0.1 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 29 0.00001 1 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 29 0.00001 10 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 50 0.00001 0.1 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 100 0.00001 0.1 1 306
echo "Line: ${LINENO} out of ${result}"
python test_model_active.py 0 0 29 1000 0.00001 0.1 1 306
