ords="2 8"
grids="12 20 28 32"
rm -f .outfile
echo "order,nxb,logerL2,logerLInf,erL2,erLInf" > .outfile
for ord in ${ords}
do
    for grid in ${grids}
    do
        ./program.opt rotate.ptl -Dord=${ord} -Dnum=${grid}
        cat .erout >> .outfile
    done
done