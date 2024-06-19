import numpy as np
import NRpyDNAcode as code
import NRpyRS as RS

import argparse
import logging
import pathlib
import subprocess
import os


def read_file(file_path):
    with open(file_path, 'r') as file:
        data = file.read().replace('\n', '')
    return data

def decode(data, error_rate, iterations):
    successful_decodes = 0
    total_bit_rate = 0

    for _ in range(iterations):
        # Introduce errors
        #error_data = code.introduce_errors(data, error_rate)
        error_data = code.createerrors(data, error_rate, 0, 0)
        # Try to decode
        #(decoded, errs_detected, errs_corrected, err_code, recoverable) = RS.rsdecode(uint8_array_length_255[, VecInt erasure_locations])"
        if is_successful:
            successful_decodes += 1
            bit_rate = len(decoded_data) / len(data)
            total_bit_rate += bit_rate
    
    if successful_decodes > 0:
        average_bit_rate = total_bit_rate / successful_decodes
    else:
        average_bit_rate = 0

    return average_bit_rate

#to be rewritten (the idea is to test all the sub error rates and then for every error rate add indels and test again)

def update_err_arr(error_array, sub_only, should_stop):
    """
    this function should update the error_array with the new values of the error rate
    """
    if sub_only:
        error_array[0] = 0.01
        sub_only = False
        #print("Substitution rate reset to 0.01")
    else:
        error_array[0] += 0.001
        if error_array[0] > 0.021:
            #print("Resetting substitution rate.")
            sub_only = True
            should_stop = True
    return sub_only, should_stop

def main():
    """
    this function has to call 100 times the test_program.py for every update of error_array or until stop condition

    """

    base_path = pathlib.Path(__file__).parent.parent.resolve()
    logger = logging.getLogger('HEDGES_benchmark_logger')
    logger.setLevel(logging.DEBUG)  # Capture all levels of messages
    console_handler = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    error_array = [0.01, 0.01, 0.01]
    should_stop = False
    sub_only = True
    first_run = True

    logger.info('Starting benchmark')
    while should_stop == False:
        py_command = ['python', str(base_path / 'python' / 'test_program.py'), '-s', str(error_array[0]), '-d',str(error_array[1]), '-i',str(error_array[2])]
        nb_success = 100
        for i in range(100) :
            process = subprocess.Popen(py_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output,error = process.communicate()
            process.wait()
            print(output)
            print(error)
            if process.returncode != 0:
                nb_success = nb_success - 1
                print("Error in subprocess")
                logger.error('Error in subprocess')
            else :
                print("Success in subprocess {}".format(i))
                logger.info('Success in subprocess {}'.format(i))
            try :
                os.rename("results/00005_560x888_94.jxl", "results/00005_560x888_94_{}.jxl".format(i))
            except :
                print("No file to remove")
        #success_rate = ((float)nb_success / 100)*100
        print("Success rate: {}% | s:{} d:{} i:{} |".format(nb_success, error_array[0], error_array[1], error_array[2]))
        logger.info('Success rate: {}% | s:{} d:{} i:{} |'.format(nb_success, error_array[0], error_array[1], error_array[2]))
        #sub_only, should_stop = update_err_arr(error_array, sub_only, should_stop)
        should_stop = True
        #update_coderate()
    logger.info('Benchmark finished')


if __name__ == '__main__':
    main()