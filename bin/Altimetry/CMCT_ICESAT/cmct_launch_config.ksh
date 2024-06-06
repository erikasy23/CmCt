## NAME: cmct_launch_config.ksh
##
## PURPOSE: Configuration for cmct_launch.ksh.  This file is *sourced* (not
## run!) by that script.  NOTE: Expects that $CMCTBIN has already been set.
##
## Configured for ggsghpcc.sgt-inc.com.
##
## AUTHOR: Jeff Guerber. SigmaSpace/GSFC 615, Mar. 2016
## 2016-04-22 JRG: Many more variables.
## 2016-04-25 JRG: BCCEMAIL. ADMINEMAIL now for error notification.
## 2016-04-26 JRG: Added Matt to BCCEMAIL and ADMINEMAIL.
## 2016-05-02 JRG: Added the new secure FTPLOCAL and FTPURL, per Craig.
## 2016-06-07 JRG: Added DLLIST.
## 2016-07-03 JRG: Changed DLLIST from ggsghpcc/dev/cmct to ggsghpcc-dev/cmct.
## 2016-07-11 JRG: Added ICESAT_CMCTNC_GRN_CFG.  Lori wants DLLIST to point
##    to ggsghpcc not ggsghpcc-dev.  Commented out FTPURL.
##
## Last SVN commit: $Id: cmct_launch_config.ksh 109 2016-07-12 07:03:38Z jguerber $

## Run directories are here
RUNDIRS=/home/cmct/RUNS

## Completed runs for retrieval.  FTPLOCAL is where we put them on the
## local file system.  DLLIST is the user's download file list.
FTPLOCAL=/WWW/secure/pub/cmct               # secure site
DLLIST=https://ggsghpcc.sgt-inc.com/cmct/download/filelist.php

## Older versions of above for reference.  FTPURL is where individual
## download packages would appear on the Web, BUT Craig says it's a
## security hole, so we're not using it anymore.
# FTPLOCAL=/home/ftp/pub/cmct               # original insecure site
# FTPLOCAL=/home/cmct/tmp-ftp/pub/cmct      # Temporary for testing
# FTPURL=https://ggsghpcc.sgt-inc.com/pub/cmct     # secure site
# FTPURL=http://ggsghpcc.sgt-inc.com/pub/cmct    # original insecure site
# DLLIST=https://ggsghpcc-dev.sgt-inc.com/cmct/download/filelist.php

DAYS2GET=10           # approx number of days available to retrieve

## BCCEMAIL will send a BCC: of the notification email to the specified
## addresses.  ADMINEMAIL are notified if anything goes wrong.
##"ES 2017-06-21: reconfigured for CmCt reinstall."

BCCEMAIL="erika.g.simon@nasa.gov"
ADMINEMAIL="erika.g.simon@nasa.gov"

## Various programs we use
JQ=/opt/csw/bin/jq      # cmd-line JSON parser
TAR=/opt/csw/bin/tar               # gnu tar
GZIP=/usr/bin/gzip

## CMCT configuration files
DATASETS_CFG=${CMCTBIN}/datasets.json
ICESAT_CISM_CFG=${CMCTBIN}/icesat_cism_config.json
ICESAT_CMCT_CFG=${CMCTBIN}/icesat_cmct_config.json    # not currently in use
ICESAT_CMCTNC_GRN_CFG=${CMCTBIN}/icesat_cmctnc_grn_config.json
CONFIGS="${DATASETS_CFG} ${ICESAT_CISM_CFG} ${ICESAT_CMCTNC_GRN_CFG}"
