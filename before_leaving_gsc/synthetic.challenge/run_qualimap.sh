for set_dir in set3 set4 set5; do
    for cell_type in normal tumour; do
	qualimap_dir=${set_dir}/qualimap/$cell_type
	mkdir -p ${qualimap_dir}
	time qualimap bamqc \
	     -bam ${set_dir}/bam/${cell_type}.bam \
	     -outdir ${qualimap_dir} \
	     -nr 10000000 \
	     -nt 20 \
	     --java-mem-size=100G
    done
done
