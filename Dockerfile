# Get repo with latest python repos
FROM python:latest

LABEL maintainer="jtuck@ncsu.edu"

WORKDIR /DORIS

COPY . .

# Install the needed tools
RUN pip3 --no-cache-dir install -r requirements.txt

# Set python path to include the code in the doris module
ENV "PYTHONPATH" "/DORIS/doris"
