# Quality analisis of run samples

Este minipipeline se ha realizado para automatizar el proceso de: **Traspaso de muestras a S3** (Por Muestras en local o en disco duro), **Cálculo de reads por muestra**, **Creación de archivo csv de muestras** (tanto grupales como individuales) y **Análisis de control de muestra por tipo de muestra**.

Todos ellos se realizaran secuencialmente usando un script de bash añadiendole los parametros necesarios para ejecutar el resto de scripts en este.

## Movimiento de Ficheros a S3


Este script mueve los archivos fastq desde una localización hasta s3 organizando todo de forma uniforme entre runs y consta de 2 posibles parametros:

- `-d` o `--date` ; Este parámetro sirve para decir la fecha del run que queremos subir a S3, es decir será desde la carpeta de ese run (local o en disco duro) y se subira al bucket de S3 `s3://raw-illumina-runs`, dentro de este se creará una carpeta con el nombre `RUN_DATE' y dentro una carpeta `FASTQ` y otra llamada `FILES`. En la carpeta `FILES` se almacenan archivos con informacion del run, mientras que en la carpeta `FASTQ`se crean unas carpetas con el formato `DATE_MUESTRA` que contrendran los archivos fastq comprimidos en gzip (tanto el *forward* como el *reverse*)

- `-w` o `--Warehouse` ; Este parámetro lo que nos permite es seleccionar la ubicacion de la carpeta el run que queremos subir. Si no se utiliza este parámetro se asume que esta en el disco duro y de allá se extraerá si por el contrario está en local hay que poner `Local` para que se coja de este lugar.

## Cálculo de reads por muestra

Este script cuenta el numero de lineas de un archivo fastq.gz y lo divide entre 4 para contabilizar el numero de reads de la muestra (se contabiliza solo los archivos *forward* ya que el numero de reads es igual en ambos archivos), para ello se utilizaran 

- `-d` o `--date` ; Este parámetro sirve para decir la fecha del run que se quiere analizar el numero de reads d sus muestras.

- `-w` o `--Warehouse` ; Este parámetro lo que nos permite es seleccionar la ubi
cacion de la carpeta el run que queremos subir. Si no se utiliza este parámetro se asume que esta en el disco duro y de allá se extraerá si por el contrario está en local hay que poner `Local` para que se coja de este lugar.

- `-p` o `--path` ; Este parámetro lo que nos genera es la ubicacion donde se guardará los archivos generados.

## Creación de archivo csv de muestras

Aqui lo que haremos es que a partir del archivo `.biom` que se ha generado en el pipeline (jenkins), se extraiga la informacion de taxonomia y de abundancia en reads para el csv global y en porcentajes para los csv individuales.




#######################################################DOCKER###################

To build Docker container is needed 
1. Storage in this folder the files with github credentials to access Run4Jenkins repository

*github_rsa*
```
-----BEGIN RSA PRIVATE KEY-----
**********************
**********************
**********************
**********************
**********************
**********************
-----END RSA PRIVATE KEY-----

```
*github_rsa.pub*
```
ssh-rsa ************************
```
2. use the command to build the image in localhost and can run a container

`docker build -t luigi/autockeck:0.10 --build-arg ssh_prv_key="$(cat ./github_rsa)" --build-arg ssh_pub_key="$(cat ./github_rsa.pub)" .`

This load the image in local computer. Then to make run the program is necessary

1. Create a file names env.list that contain the credentials for different conexions
```
# AWS Credentials
AWS_ACCESS_KEY_ID=********************
AWS_SECRET_ACCESS_KEY=********************
# BS Credentials
ClientKey=********************
ClientSecret=********************
AccessToken=********************
# DB Credentials
DB_USER=********************
DB_PASSWORD=********************
DB_HOST=********************
DB_PORT=********************
DB_NAME=********************

```
2. Run the image of container with parameters of env.list and RunDate

`docker run --name Run20181113 --env-file env.list -e RunDate=20181113 luigi/autockeck:0.10`

--name is the container name that will be used to have a pattern to copy from container the result folder

The last step is take the results to localhost for this the command will be.

`docker cp Run20181113:/Results/20181113/ ./Results/`
