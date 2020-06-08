#!/bin/bash

hmodel=SIBYLL
inbase=/data/sim/IceCube/2011/generated/dcorsika/length_1600_radius_800_v6960-5comp_sibyll_5-component/emin_600_emax_1e11_dslope_0_pgam_E2.0E2.0E2.0E2.0E2.0
outbase=/data/uwa/jvansanten/projects/2012/muongun/corsika/$hmodel

for atmod in 12 11 13 14; do
  for flux in 5Comp CascadeOptimized Hoerandel5 GaisserH3a GaisserH4a; do
      outdir=$outbase/$flux/atmod_$atmod
      if [[ ! -d $outdir ]]; then
        mkdir -p $outdir
      fi
      
      tags=""
      outfiles=""
      # fill file into a histogram
      # redirected loop: while; do .. done < <(command) is bash-only syntax
      # for running a loop without pipes. the piped loop would run in a
      # subshell and be unable to modify variables in the parent loop
      i=0
      while read infiles; do
        basename=corsika_2500000.$(printf %.6d $i)
        outfile=$outdir/$basename.hdf5
        tag=histogram.$basename.$hmodel.atmod_$atmod.$flux
        echo JOB $tag fill_histograms.sub
        echo VARS $tag args=\"$infiles $outfile --flux=$flux --detcfg=1 --mindepth=1 --maxdepth=2.8 --steps=19\"
        tags="$tags $tag"
        outfiles="$outfiles $outfile"
        i=$(($i+`echo $infiles | wc -w`))
      done < <(find $inbase/atmod_$atmod -iname 'corsika_2500000_*.gz' | sort | xargs -L100 echo)
      # now, merge them all together
      tag=histogram_merge.$hmodel.atmod_$atmod.$flux
      echo JOB $tag histadd.sub
      echo VARS $tag arguments=\"--overwrite --norm=$i $outfiles $outdir.hdf5\"
      echo PARENT $tags CHILD $tag
      break
  done
  break
done
