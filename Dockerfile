# Stage 1: Build the C++ code
FROM python:2.7 AS build

WORKDIR /app

COPY . /app

RUN apt-get update && apt-get install -y \
    build-essential \
    python-dev

RUN pip install numpy

RUN pip install pathlib

RUN make

RUN cp ./cpp/NRpyRS.so ./cpp/NRpyDNAcode.so $(python -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())")

# Keep the container running
CMD ["tail", "-f", "/dev/null"]