#!/bin/bash -f
# this is a script for processing LIFETRAC output file .ltr
# to produce the following step-by-step data:
#     - normalized beam intensity: 'intensity.txt'
#     - normalized rate of losses: 'lossrate.txt'
#     - specific luminosity (L/L0*N/N0): 'luminosity.txt'
#     - horizontal and vertical emittances: 'emit.txt'
#     - longitudinal beam sigma: 'sigm.txt'
#
if [ $# != 3 ]; then
    if [ $# == 1 ]; then
         echo "Interaction Points for which the luminosity is monitored are:"
         grep -A 3 Luminosity $1 | head -3
         exit 0;
    fi
    echo "This is a script for processing LIFETRAC output"
    echo " "
    echo "invocation:"
    echo "     $0 ltr_file Watch_IP py_dir"
    echo " "
    echo "where Watch_IP is the name of the IP to print luminosity, e.g. \"IP_1077\"."
    echo "and py_dir is the path to the python LifeDesk module"
    echo " "
    echo "To find out which IPs are in your ltr file run: "
    echo "$0 ltr_file"
    exit 0;
fi

PYDIR=$3

grep 'Mean_values' $1 | cut -f2 -d'=' > inti.txt
${PYDIR}/normint < inti.txt > intensity.txt
${PYDIR}/lossrate2 < intensity.txt > lossrate.txt
rm inti.txt

grep ": $2" $1 | cut -f3 -d' ' > lumi.txt
${PYDIR}/speclumi > luminosity.txt
rm lumi.txt

${PYDIR}/emitgfit $1 > emit.txt
grep '|sigm|' $1 | cut -f10 -d' ' > sigm.txt
