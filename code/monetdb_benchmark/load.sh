#!/bin/sh
set -x
if [ "$#" -ne 4 ]; then
    echo "Usage: load.sh DATA_PATH NGENES NPATIENS DATABASE"
	exit -1
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DPATH=$1
NGENES=$2
NPATIENTS=$3
DATABASE=$4
GEO="${DPATH}/GEO-${NGENES}-${NPATIENTS}.txt"
GO="${DPATH}/GO-${NGENES}-${NPATIENTS}.txt"
GENES="${DPATH}/GeneMetaData-${NGENES}-${NPATIENTS}.txt"
PATIENTS="${DPATH}/PatientMetaData-${NGENES}-${NPATIENTS}.txt"


mclient $DATABASE $DIR/make.sql

# Load the data
mclient $DATABASE -s "COPY OFFSET 2 INTO geo FROM '$GEO' USING DELIMITERS ', ' LOCKED;"
mclient $DATABASE -s "COPY OFFSET 2 INTO go_matrix FROM '$GO' USING DELIMITERS ', ' LOCKED;"
mclient $DATABASE -s "COPY OFFSET 2 INTO genes FROM '$GENES' USING DELIMITERS ', ' LOCKED;"
mclient $DATABASE -s "COPY OFFSET 2 INTO patients FROM '$PATIENTS' USING DELIMITERS ', ' LOCKED;"

