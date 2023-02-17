for f in pacific_*grid_T*nc;
do
    echo "Processing $f"

    if [ -f "compressed_$f" ]; then
        echo "$p processed. skip"
        continue
    fi
    ncks -L 9 -v thetao,e3t $f compressed_$f
done

for f in pacific_*speed_U*nc;
do
    if [ -f "compressed_$f" ]; then
        echo "$p processed. skip"
        continue
    fi
    echo "Processing $f"
    ncks -L 9 -v uo $f compressed_$f
done

for f in pacific_*speed_V*nc;
do
    if [ -f "compressed_$f" ]; then
        echo "$p processed. skip"
        continue
    fi
    echo "Processing $f"
    ncks -L 9 -v vo $f compressed_$f
done

for f in pacific_*add_T*nc;
do
    if [ -f "compressed_$f" ]; then
        echo "$p processed. skip"
        continue
    fi
    echo "Processing $f"
    ncks -L 9 -v NCHL,DCHL $f compressed_$f
done

for f in pacific_*ptrc_T*nc;
do
    if [ -f "compressed_$f" ]; then
        echo "$p processed. skip"
        continue
    fi
    echo "Processing $f"
    ncks -L 9 -v ZOO2,ZOO,GOC,PHY2 $f compressed_$f
done
