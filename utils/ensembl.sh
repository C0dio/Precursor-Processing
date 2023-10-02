#!/bin/bash

#
# Default Variables
#
SEQUENCE_ID="${SEQUENCE_ID:-ENSG00000139618}"
SPECIES="${SPECIES:-""}"
TYPE="${TYPE:-""}"
GFF3_OUTPUT="${GFF3_OUTPUT:-inputs/input.gff3}"
FASTA_OUTPUT="${FASTA_OUTPUT:-inputs/input.fa}"

set -euo pipefail

#
# Validate the Sequence ID
#
LOOKUP=$(curl "https://rest.ensembl.org/lookup/id/${SEQUENCE_ID}?expand=0;format=condensed" -H "Content-type:application/json")
if [[ "$LOOKUP" == *"error"* ]]; then
    echo "ERROR: Invalid Sequence Id: $SEQUENCE_ID"
    exit 1
fi

#
# Get the species
#
if [ -z "$SPECIES" ]; then
    SPECIES=$(echo $LOOKUP | grep -oP '(?<="species":)".+?"' | tr -d '"')
fi

#
# Get the sequence type
#
if [ -z "$TYPE" ]; then
    TYPE=$(echo $LOOKUP | grep -oP '(?<="object_type":)".+?"' | tr -d '"')
fi

#
# Get the GFF3 data
#
curl -L "https://www.ensembl.org/${SPECIES}/Export/Output/${TYPE}?g=${SEQUENCE_ID};output=gff3;param=gene;param=transcript;param=intron;_format=Text" -o $GFF3_OUTPUT

#
# Validate file
#
# if HTML was returned then the endpoint failed to format to text
#
TEMP=$(head $GFF3_OUTPUT -n 1)
if [ "$TEMP" == "<!DOCTYPE html>" ]; then
    echo "ERROR: File Validation Failed - Species was Invalid"
    echo "ERROR: Sequence Id: $SEQUENCE_ID"
    echo "ERROR: Object Type: $TYPE"
    echo "ERROR: Species: $SPECIES"
    exit 1
fi

#
# Get the FASTA data
# CL TODO: Uses the same URL as GFF3 but different format (DRY Issue)
#
curl -L "https://www.ensembl.org/${SPECIES}/Export/Output/${TYPE}?g=${SEQUENCE_ID};output=fasta;param=gene;param=transcript;param=intron;_format=Text" -o inputs/input.fa