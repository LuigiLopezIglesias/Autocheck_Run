# Use an official Python runtime as a parent image
FROM luigi/autockeck:0.05 

## Instal repositories
#RUN apt-get update \
#    && apt-get upgrade -y \
#    && apt-get install -y \
#    apt-utils \
#    build-essential \
#    ca-certificates \
#    gcc \
#    git \
#    libpq-dev \
#    make \
#    python-pip \
#    python2.7 \
#    python2.7-dev \
#    ssh \
#    libc6 \
#    libgcc1 \
#    libstdc++6 \
#    python-cogent \
#    python-dateutil \
#    biom-format-tools \ 
#    && apt-get autoremove \
#    && apt-get clean \
#    && pip install --upgrade pip 
#
#RUN pip install python-dateutil==2.7.5
 
COPY requirements.txt Analysis.py / 
#RUN pip install -r requirements.txt

COPY PyScripts/* /PyScripts/
COPY github_rsa /.shh/github_rsa

# Run app.py when the container launches
CMD ["python", "/Analysis.py", "-d", "20181113"]
