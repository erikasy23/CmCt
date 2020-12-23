#!/bin/sh

SUBPROGRAM=CMCT_ICESAT

# needed by the CmCt program
CMCTBIN=/opt/cmct/${SUBPROGRAM}

# For now we will be sloppy and use /opt/cmct/CMCT_ICESAT as working/compiling/etc.
# Once you run the actual docker container you go in and do the following:
rsync -avz /scratch/${SUBPROGRAM} /opt/cmct && \
cd ${CMCTBIN} && \
make all

