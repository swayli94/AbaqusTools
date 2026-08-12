#!/bin/bash
#
# Run a wingbox case in the current folder.
#
#   ./run.sh                             # parameters-simple-wingbox.json
#   ./run.sh parameters-inner-wingbox.json    # any other parameter file
#
# The working folder must contain, next to this script:
#
#   AbaqusTools/        the package folder (from the repository root)
#   *.py                every python file of `wingbox/`, i.e., run.py,
#                       wingbox_model.py, params.py, postprocess_failure.py,
#                       extract_results.py, lofting_part.py, rib_part.py,
#                       geometry.py, layup.py, utils.py
#   BAC-NLF.dat         the airfoil referenced by `pGeo.sections[*].airfoil`
#   <parameters>.json   the parameter file passed to this script
#   *.f90               only for "failure_model": "LaRC05", the user
#                       subroutines of `LaRC05/` (uvarm.f90 and its includes)
#
PARAMS="${1:-parameters-simple-wingbox.json}"

Here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${Here}"

if [ ! -f "${PARAMS}" ]; then
    echo "Parameter file not found: ${PARAMS}"
    exit 1
fi

# Clean the results of the previous run.  Parameter files, python files, the
# airfoil and the user subroutines are kept.
rm -f *.pyc *.rec *.exception fort* abaqus.rpy* *.cae *.jnl *.dbg *.rpy
rm -f abaqus_acis.log abq.app_cache abaqus.rpt
rm -f Job_*.odb Job_*.dat Job_*.msg Job_*.sta Job_*.log Job_*.prt Job_*.inp
rm -f Job_*.com Job_*.env Job_*.sim Job_*.stt Job_*.lck Job_*.odb_f
rm -f Job_*_failure_envelope.npz Job_*_failure_summary.json Job_*_mass.json
rm -f Job_*.SMABulk Job_*.dmp.lnz.*
rm -f larc05_fi_track_*.txt

echo "wingbox: ${PARAMS}"

python -u run.py --params "${PARAMS}"
