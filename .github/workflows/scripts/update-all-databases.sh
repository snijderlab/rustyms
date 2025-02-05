#!/usr/bin/env bash
set -euo pipefail


function help {
    echo "Usage: generate-all-databases.sh"
    echo ""
    echo "Download the required databases and build the "
    echo "required binary representations of the ontologies."
    echo ""
    echo "Options:"
    echo "  -h, --help     Display this help and exit"
    exit 1
}

# Download all of the relevant databases, decompressing, filtering, and moving
# them as needed.
function download-dbs {
    echo "Downloading IMGT..."
    mkdir -p rustyms-generate-imgt/data
    curl https://www.imgt.org/download/LIGM-DB/imgt.dat.Z \
        | gunzip -c > rustyms-generate-imgt/data/imgt.dat


    echo "Downloading databases..."
    db_data="rustyms-generate-databases/data"
    mkdir -p ${db_data}
    curl https://raw.githubusercontent.com/HUPO-PSI/psi-mod-CV/refs/heads/master/PSI-MOD-newstyle.obo \
        > ${db_data}/PSI-MOD-newstyle.obo
    curl http://www.unimod.org/obo/unimod.obo > ${db_data}/unimod.obo
    curl ftp://ftp.proteininformationresource.org/pir_databases/other_databases/resid/RESIDUES.XML \
        > ${db_data}/RESID-RESIDUES.XML
    curl https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/cv/XLMOD.obo \
        > ${db_data}/XLMOD.obo
    curl -L http://purl.obolibrary.org/obo/gno.obo \
        | sed '/(property_value: GNO:00000(022|023|041|042|101|102) .*$\n)|(def: .*$\n)/d' \
        | gzip -c \
        > ${db_data}/GNOme.obo.gz
    curl https://glycosmos.org/download/glycosmos_glycans_list.csv \
        | cut -f1,2 -d',' \
        | gzip -c > ${db_data}/glycosmos_glycans_list_filtered.csv.gz
}

# Serialize the databases to binary blobs to build into rustyms.
function build-binaries {
    echo "Serializing IMGT..."
    cargo run --bin rustyms-generate-imgt
    echo "Serializing the other databases..."
    cargo run --bin rustyms-generate-databases
}

function main {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -h|--help)
                help
                ;;
            *)
                echo "Unknown argument: $1"
                help
                ;;
        esac
    done

    download-dbs
    build-binaries
}

main "$@"
