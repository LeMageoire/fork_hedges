import numpy as np
import NRpyDNAcode as code
import NRpyRS as RS

import argparse
import logging
import pathlib


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

def main():
    base_path = pathlib.Path(__file__).parent.parent.resolve()
    data = read_file(str(base_path / "data" / "D"))
    average_bit_rate = decode(data, 0.02, 100)
    print("Average bit rate: {}".format(average_bit_rate))

if __name__ == '__main__':
    main()