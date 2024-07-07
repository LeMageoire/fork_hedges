import numpy as np
import NRpyDNAcode as code
import NRpyRS as RS

import os
import argparse
import logging
import pathlib
import subprocess


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
    args = argparse.ArgumentParser(description='Benchmark for HEDGES, a DNA storage encoder/decoder')
    args.add_argument('-s', '--substitution', type=float, help='Substitution error rate')
    args.add_argument('-d', '--deletion', type=float, help='Deletion error rate')
    args.add_argument('-i', '--insertion', type=float, help='Insertion error rate')
    args.add_argument('-r', type=int, default=100, help='Number of iterations')
    args.add_argument('-t', '--target', type =str, help='input file path')
    args.add_argument('-o', '--output', type =str, default = "./results/default", help='folder_path + file path')
    args = args.parse_args()

    # Check if the input file can be opened
    try:
        with open(args.target, 'r') as f:
            pass  # Just checking if the file can be opened
    except IOError:
        logger.error(f"Unable to open the input file: {args.target}")
        return
    
    # Check if the output file can be opened (create the folder if it doesn't exist)
    
    output_dir = pathlib.Path(args.output).parent
    try:
        os.makedirs(output_dir, exist_ok=True)
    except IOError:
        logger.error(f"Unable to open the output folder: {args.output}")
        return

    base_path = pathlib.Path(__file__).parent.parent.resolve()
    logger = logging.getLogger('HEDGES_benchmark_logger')
    logger.setLevel(logging.DEBUG)  # Capture all levels of messages
    console_handler = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    if args.substitution is None or args.deletion is None or args.insertion is None:
        error_array = [0.017, 0.017, 0.017] # Manual Configuration
    else :
        error_array = [args.substitution, args.deletion, args.insertion]
    should_stop = False
    sub_only = True
    first_run = True

    logger.info('Starting benchmark')
    while should_stop == False:
        py_command = ['python', str(base_path / 'python' / 'test_program.py'), '-s', str(error_array[0]), '-d',str(error_array[1]), '-i',str(error_array[2], '-r', str(args.r), '-t', str(args.target), '-o', str(args.output)]
        nb_success = 100
        output_path = Path(args.output)
        base_name = output_path.stem
        extension = output_path.suffix
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
                new_file_name = "{}{}_{}{}".format(output_path, base_name, i, extension)
                os.rename(args.output, new_file_name)
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