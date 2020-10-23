#!/bin/ksh -x
##
## NAME:cmct_launch.ksh
##
## PURPOSE: Launch script for CMCT.  Gets information about a new
## submission, creates a run directory, populates it with the downloaded
## model file and .json configuration files, and launches cmct_main.
##
## The run directory will have 3 subdirectories, MODELS, WORK, and OUT.
## Model files go in MODELS.  The work is done in WORK.  Output is packaged
## into OUT before being copied to the retrieval directory.
##
## REQUIRES:  The "jq" command-line JSON parser.
##
## PARAMETERS:
##    ${1}:  Path to JSON run configuration file
##
## TBD: Error checking is minimal!!!
##
## AUTHOR: Jeff Guerber, SigmaSpace/GSFC615, Mar. 2016
## 2016-04-21 JRG: Latest.
## 2016-04-22 JRG: Tar and gzip the results, copy to "ftp" area, send user
##    email.  Bug fixes.
## 2016-04-25 JRG: Added error trapping.
## 2016-04-26 JRG: chmod g+w on created directories.  Bug fixes.
## 2016-04-27 JRG: Verbose flag -v on mailx calls. Fix BCC.
## 2016-04-27 JRG: -r cmct@sgt-inc.com on mailx calls.
## 2016-05-02 JRG: Added security code from Craig.
## 2016-06-02 CJ?: Updated .htaccess content.
## 2016-06-07 JRG: Copy result files to a staging dir and tar that.  Added
##    the new results directory list URL to the user email. (Direct
##    retrieval seems to still work too.)
## 2016-06-07 JRG: Lori doesn't want direct link ($URL) in the email.
##    Added file name and run title to email.
## 2016-06-10 JRG: At Tom's request, put the direct link back in the email.
## 2016-06-10 JRG: Turns out the direct link is vulnerable to SQL
##    injection attacks!  So it's back out, again.
## 2016-07-11 JRG: Minor change in the email text.
##
## Last SVN commit: $Id: cmct_launch.ksh 108 2016-07-12 06:59:10Z jguerber $
##"ES 2017-06-21: reconfigured for CmCt reinstall."

## Source configuration, from same directory as this script.
CMCTBIN=$(dirname ${0})
. ${CMCTBIN}/cmct_launch_config.ksh

RUNCFG=${1}

## Function traperror handles any errors that occur, so the script doesn't
## just try to continue.  Sends email to $ADMINEMAIL and if possible to the
## user, then exits with status 1.
traperror()
{
    status=$?
    echo "Error trap taken, status=${status}.  Sending emails."

    mailx -v -s "CMCT Error occurred in cmct_launch.ksh" -r cmct@sgt-inc.com ${ADMINEMAIL} <<EOF
An error has occurred running cmct_launch.ksh on ggsghpcc.sgt-inc.com.

Run config file = ${RUNCFG}
Status = ${status}

Please check the launch log.
EOF

    if [[ -n ${USEREMAIL} ]]
    then
	# This should probably include contact information
	mailx -v -s "CMCT error runid=${RUNID}" -r cmct@sgt-inc.com ${USEREMAIL} <<EOF
${REALNAME}:

    We're sorry, an error occurred during processing of your CMCT job, runid = ${RUNID}.  We've notified the CMCT staff to take a look.

    This is an automated email.  Please don't reply to it directly.
EOF
    fi

    exit 1
}

##
## Main part of the script
##

trap 'traperror' ERR

RUNCFG_BASE=$(basename ${RUNCFG})

##
## Parse what we need from the run config file
##
RUNID="$(      ${JQ} -r '.cmct_run_config.run.runid' < ${RUNCFG} )"
UPLDIR="$(     ${JQ} -r '.cmct_run_config.run.upload_dir' < ${RUNCFG} )"
LOGINNAME="$(  ${JQ} -r '.cmct_run_config.run.loginname' < ${RUNCFG} )"
REALNAME="$(   ${JQ} -r '.cmct_run_config.run.actualusername' < ${RUNCFG} )"
USEREMAIL="$(  ${JQ} -r '.cmct_run_config.run.email' < ${RUNCFG} )"
RUNTITLE="$(   ${JQ} -r '.cmct_run_config.run.user_run_title' < ${RUNCFG} )"
MODELFILES="$( ${JQ} -r '.cmct_run_config.comparisons[].model.filename' < ${RUNCFG} )"

##
## Make the run dirs.
## Model files go in subdir MODELS
##
RUNDIR=${RUNDIRS}/${RUNID}
WORKDIR=${RUNDIR}/WORK
MODELDIR=${RUNDIR}/MODELS
OUTDIR=${RUNDIR}/OUT

mkdir ${RUNDIR}
mkdir ${WORKDIR} ${MODELDIR} ${OUTDIR}

# Add g+w so maybe cmct staff can clean up.  (But with sticky bit set,
# maybe we can't anyway.)
chmod g+w ${RUNDIR} ${WORKDIR} ${MODELDIR} ${OUTDIR}

##
## Change to work dir.  Create the control file.  Populate the work dir.
##
cd ${WORKDIR}
CTLFILE=${WORKDIR}/cmct_main_control_${RUNID}.json  # Name of the control file
cat  > ${CTLFILE}  <<EOF
{
    "COMMENT" : [
	"CMCT Main Program control file.",
        "Written by cmct_launch.ksh",
        "$(date)"
    ],

    "cmct_main" : {
	"working_dir" : "${WORKDIR}",
        "models_dir"  : "${MODELDIR}",
	"cmct_run_config" : "${RUNCFG_BASE}",
	"verbosity" : 0
    }
}
EOF

cp -p ${RUNCFG} ${WORKDIR}
cp -p ${CONFIGS} ${WORKDIR}

# Copy models from UPLDIR to MODELDIR
for mf in ${MODELFILES}
do
    cp -p ${UPLDIR}/${mf} ${MODELDIR}
done

##
## Launch the program.  Put both stdout and stderr in the log file.  Don't
## trap non-zero exit running cmct_main, we want to report even errors to
## the user, but turn trapping back on afterwards.
##
trap - ERR
LOGFILE=${WORKDIR}/CMCT_${RUNID}.cmct_main.log
${CMCTBIN}/cmct_main ${CTLFILE} > ${LOGFILE} 2>&1
trap 'traperror' ERR

##
## Tar up and gzip the results into $OUTDIR.
## Copy into a staging directory, then tar that, so that we get a top-level
## directory in the tar file.
## We may eventually be more selective about what we send.
##
STAGEDIR=CMCT_${RUNID}
mkdir ${OUTDIR}/${STAGEDIR}
cp -pr * ${OUTDIR}/${STAGEDIR}

TARFILE=${OUTDIR}/CMCT_${RUNID}.tar
${TAR} cvf ${TARFILE} -C ${OUTDIR} ${STAGEDIR}
${GZIP} --best ${TARFILE}
TGZFILE=${TARFILE}.gz             # assumes std gzip name convention

rm -rf ${OUTDIR}/${STAGEDIR}

##
## Copy it to the "ftp" area for user retrieval (even though it's really
## done by html not ftp, but Craig the sysadmin still calls it that).
##
USERFTP=${FTPLOCAL}/${LOGINNAME}

# Check that user's directory exists, if not (probably a new user) create it.
if [[ ! -d ${USERFTP} ]]
then
    mkdir ${USERFTP}
fi
chmod g+rw ${USERFTP}          # Let group cmct have access (for now anyway)

# Security code from Craig.  Requires user to log in, and doesn't let them
# see anyone else's files.
if [ ! -f $USERFTP/.htaccess ]
then
   echo "
AuthType Basic
AuthName \"Files for $LOGINNAME\"
AuthBasicProvider dbd

require user $LOGINNAME" > $USERFTP/.htaccess
fi

cp -p ${TGZFILE} ${USERFTP}

# $URL is the direct link to the results file.
# I guess its vestigal now.  Used to include it in the email, but Craig
# says it's vulnerable to SQL injection attacks.
# URL=${FTPURL}/${LOGINNAME}/$(basename ${TGZFILE})

##
## Email the user that the files are ready
##
## $BCCEMAIL must be in quotes here, else multiple BCC addresses may be
## treated as regular ones.  Seems to be OK if it's blank.
mailx -v -s "Your CMCT results are ready: ${RUNID}" -r cmct@sgt-inc.com -b "${BCCEMAIL}" ${USEREMAIL} <<EOF
${REALNAME}:

   Thank you for using the ISMIP6/NASA GSFC Cryosphere Model Comparison Tool.

   The results of your recent run are ready!  For at least the next ${DAYS2GET} days, you can retrieve them from your personal download directory at this URL:

${DLLIST}

   You will need to log in with the userid and password you use on the CMCT submission site.

   Run title:  ${RUNTITLE}
   Run's unique runid:  ${RUNID}
   File name of run results package:  $(basename ${TGZFILE})

   If you have any problems or questions or comments, please contact:

cmct@sgt-inc.com                              (general contact address)
Sophie Nowicki <sophie.nowicki@nasa.gov>      (project leader)
Erika Simon <erika.g.simon@nasa.gov>          (software and comparison help)
Lori Tyahla <LTyahla@sgt-inc.com>             (web site help only)

   This is an automated email.  Please don't reply to it directly.
EOF

echo "Email notification sent to ${USEREMAIL} for runid ${RUNID}"
echo
times

exit 0
