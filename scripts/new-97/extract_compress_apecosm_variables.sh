for p in pacific_ORCA1*nc;
do
    fout="compressed_$p"
    if [ -f "$fout" ]; then
        echo "$p processed. skip"
        continue
    fi

    echo $p
    ncks -L 9 $p $fout
done
