d="run/s03"; rm -rf $d/{cms,std/log}/*; mkdir -p $d
for i in {2001..3000}; do
    for m in {10..50..10}; do
	for f in {05,15,25}; do
            echo -n "Rscript -e 'source(\"ini.R\"); "
	    echo    "r <- sim(N=1500, M=$m, NAF=.$f, tag=\"F$f\", seed=$i, out=\"{n:05d}.mix\")'"
	done
    done
done | hpcwp - -d$d -q45 -m4 -p1 --wtm 4 --ln '*.rds' --cp "*.R" --cp "R" --tag ${d##*/}
d="run/s01"; Rscript -e "source('rpt.R'); plt.sim('$d')"


d="run/mc7"; rm -rf $d/{cms,std/log}/*; mkdir -p $d
for i in {1001..2000}; do
    # for m in {10..50..10}; do
    for f in {05,15,25}; do
        echo -n "Rscript -e 'source(\"ini.R\"); "
	echo    "mic(N=1500, M=70, NAF=.$f, tag=\"F$f\", seed=$i, out=\"{n:04X}.mic\")'"
    done
    # done
done | hpcwp - -d$d -q6 -m4 -p1 --wtm 2 --ln '*.rds' --cp "*.R" --cp "R" --tag ${d##*/}

