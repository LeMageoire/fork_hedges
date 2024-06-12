#!/bin/sh

#this function move the 2 .so files into the ./venv/lib/python3.9/site-packages


echo " .so move"
mv cpp/NRpyDNAcode.so ./venv/lib/python3.9/site-packages
mv cpp/NRpyRS.so ./venv/lib/python3.9/site-packages
echo ".so done moving"