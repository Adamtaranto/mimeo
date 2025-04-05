#!/bin/sh

# Description
echo "Screen genome A (target) for features which are high copy in genome B (query)." >&2

# Usage example:
##./Interspecies_Repeat_Scan -z /usr/local/bin/lastz -t target_split -q query_split -i 60 -p query2target_60 -l 100 -d outdir -c 5 -L chromlens.txt -r

# Reset OPT index to 1
OPTIND=1

# Initialize our own variables:
minIdt="60"
minLen="100"
outdir="outdir"
prefix="genome_alignment"
coverage="5"
skipSelf=0
LZ="/usr/local/bin/lastz"
chromLens=0
recycleTAB=0

while getopts "d:p:t:q:i:l:z:c:L:rs" opt; do #Options followed by ":" expect an arguement.
	case "$opt" in
	z)	LZ=$OPTARG
		if [ ! -f "$LZ" ]; then
			echo "LASTZ does not exist at path: $LZ" >&2
			exit 1
		fi
		echo "Path to LASTZ app set as: $OPTARG" >&2
		;;
	L)	chromLens=$OPTARG
		if [ ! -f "$chromLens" ]; then
			echo "Seq lengths file does not exist at path: $chromLens" >&2
			echo "Create length file using samtools faidx:" >&2
			printf %"s\n" "samtools faidx Target_Genome.fa" "cut -f1,2 Target_Genome.fa.fai > Target_Scaffold.sizes"
			exit 1
		fi
		;;
	i)	minIdt=$OPTARG
		echo "Setting min identity threshold at $OPTARG %" >&2
		;;
	l)	minLen=$OPTARG
		echo "Setting min hit length threshold to $OPTARG" >&2
		;;
	d)	outdir=$OPTARG
		echo "Setting output directory to $OPTARG" >&2
		;;
	p)	prefix=$OPTARG
		;;
	t)	tdata=$OPTARG
		if [ ! -d "$tdata" ]; then
			echo "Target sequence directory $tdata does not exist." >&2
			exit 1
		fi
		;;
	q)	qdata=$OPTARG
		if [ ! -d "$qdata" ]; then
			echo "Query sequence directory $qdata does not exist." >&2
			exit 1
		fi
		;;
	s)	skipSelf=1
		echo "Ignoring any sequence pairs with same name." >&2
		;;
	r)	recycleTAB=1
		echo "Use existing alignment data if found." >&2
		;;
	c)	coverage=$OPTARG
		echo "Coverage threshold set as: $OPTARG" >&2
		;;
	\?)	echo "Invalid option: -$OPTARG" >&2
		exit 1
		;;
	:)
		echo "Option -$OPTARG requires an argument." >&2
		exit 1
		;;
	esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

echo "Leftover arguements:" "$@" >&2

# Housekeeping
if [ "$chromLens" == 0 ]; then
	echo "Create length file using samtools faidx:" >&2
	printf %"s\n" "samtools faidx Target_Genome.fa" "cut -f1,2 Target_Genome.fa.fai > Target_Scaffold.sizes"
	exit 1
fi

if [ ! -d "$outdir" ]; then
	echo "Creating output directory $outdir" >&2
	# Make output directory if does not exist
	mkdir "$outdir"
fi

# Make empty output files
out_gff_raw=$(echo $outdir"/"$prefix"_raw_alignments.gff3")
out_gff_cov=$(echo $outdir"/"$prefix"_covFilter_alignments.gff3")
out_tab=$(echo $outdir"/"$prefix"_concat.tab")
temp_bed=$(echo $out_tab".temp.bed")
temp_bed_sorted=$(echo $out_tab".temp.sorted.bed")

# Scrub old results
if [ -f "$out_tab" ] && [ $recycleTAB == 0 ]; then
	echo "Removing old output file: $out_tab" >&2
	rm "$out_tab"
	# Initialise output files
	echo $'#name1\tstrand1\tstart1\tend1\tname2\tstrand2\tstart2+\tend2+\tscore\tidentity' > $out_tab
elif [ ! -f "$out_tab" ]; then
	# Initialise output files
	echo $'#name1\tstrand1\tstart1\tend1\tname2\tstrand2\tstart2+\tend2+\tscore\tidentity' > $out_tab
	recycleTAB=0
else
	echo "Use existing alignment file: $out_tab" >&2
fi

# Scrub old raw hits GFF if not getting used
if [ -f "$out_gff_raw" ] && [ $recycleTAB != 1 ]; then
	echo "Removing old output file: $out_gff_raw" >&2
	rm "$out_gff_raw"
fi

# Scrub old results GFF
if [ -f "$out_gff_cov" ]; then
	echo "Removing old output file: $out_gff_cov" >&2
	rm "$out_gff_cov"
fi

# Initialise coverage GFF file
echo $'##gff-version 3\n#seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes' > $out_gff_cov

if [ $recycleTAB != 1 ]; then
	# Initialise raw hit GFF file
	echo $'##gff-version 3\n#seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes' > $out_gff_raw
	# Run alignments
	for target in $(find $tdata -type f \( -name "*.fa" -or -name "*.fasta" -or -name "*.fsa" \) | sort)
	do
		for query in $(find $qdata -type f \( -name "*.fa" -or -name "*.fasta" -or -name "*.fsa" \) | sort)
		do
		t_file=$(basename $target)
		t_name="${t_file%.*}"
		q_file=$(basename $query)
		q_name="${q_file%.*}"
		outfile=$(echo $outdir"/"$t_name".vs."$q_name".tab")
		#feature=$(echo $q_name"_ident")
		if [ $t_name == $q_name ] && [ $skipSelf != 0 ]; then
			echo "Skip comparison of $t_name and $q_name. Same sequence." >&2
			continue
		else
			$LZ $target $query \
			--gfextend \
			--chain \
			--gapped \
			--step=1 \
			--strand=both \
			--hspthresh=3000 \
			--entropy \
			--output=$outfile \
			--format=general:name1,strand1,start1,end1,length1,name2,strand2,start2+,end2+,length2,score,identity \
			--verbosity=1 \
			--markend
			#Scrub % symbols
			sed -i '' -e 's/%//g' $outfile
			## Filter Inter_Chrome targets to min len $minLen [100], min identity $minIdt [90]
			## New fields = name1,strand1,start1,end1,name2,strand2,start2+,end2+,score,identity
			## Sort filtered bed file by chrom, start, stop
			echo "Alignment finished, writing hits: $out_tab" >&2
			awk '!/^#/ { print; }' $outfile | awk -v minLen="$minLen" '0+$5 >= minLen {print ;}' | awk -v OFS='\t' -v minIdt="$minIdt" '0+$13 >= minIdt {print $1,$2,$3,$4,$6,$7,$8,$9,$11,$13;}' | sed 's/ //g' | sort -k 1,1 -k 3n,4n >> $out_tab
			## Create GFF3 file for merged filtered hits
			echo "Writing filtered hits to gff3: $out_gff_raw" >&2
			#awk '!/^#/ {print ;}' $out_tab | awk -v OFS='\t' -v q_name="$q_name" 'BEGIN{i=0}{i++;}{j=sprintf("%09d",i)}{print $1,"LASTZ","'"$feature"'",$3,$4,$9,$2,".","ID=LZ_hit_"q_name"_"j";Idt="$10";Target="$5"_"$6"_"$7"_"$8 ;}' >> $out_gff_raw
			awk '!/^#/ {print ;}' $out_tab | awk -v OFS='\t' -v q_name="$q_name" 'BEGIN{i=0}{i++;}{j=sprintf("%09d",i)}{print $1,"LASTZ","Raw_Alignment",$3,$4,$9,$2,".","ID=LZ_hit_"q_name"_"j";Idt="$10";Target="$5"_"$6"_"$7"_"$8 ;}' >> $out_gff_raw
			rm $outfile
		fi
		done
	done
fi

## Port filtered hits to bed format: "name1,start1,end1"
echo "Converting filtered hits to bed" >&2
awk '!/^#/ {print $1,"\t",$3,"\t",$4;}' $out_tab > $temp_bed
## Remove spaces
sed -i '' -e 's/ //g' $temp_bed
## Sort filtered bed file by chrom, start, stop
echo "Sorting bed file" >&2
sort -k 1,1 -k 2n,3n $temp_bed > $temp_bed_sorted
## Generate non-zero coverage scores for target genome regions, filter for min coverage of x
# Note: Recycling $temp_bed
echo "Calculate target genome coverage and filter for min $coverage X coverage. Checking for 'bedtools genomecov'..." >&2
bedtools genomecov -bg -i $temp_bed_sorted -g $chromLens | awk -v OFS='\t' -v cov="$coverage" '0+$4 >= cov {print ;}' > $temp_bed
## Re-sort
sort -k 1,1 -k 2n,3n $temp_bed > $temp_bed_sorted
## Merge end-to-end annotations
echo "Merge overlapping annotations." >&2
bedtools merge -i $temp_bed_sorted > $temp_bed
## Filter orphan fragments < xx && Create GFF3 file
echo "Writing filtered hits to gff3: $out_gff_cov" >&2
awk -v minLen="$minLen" '{ if($3 - $2 >= minLen) print ;}' $temp_bed | awk -v OFS='\t' 'BEGIN{i=0}{i++;}{j= sprintf("%05d", i)}{print $1,"LASTZ","'"$prefix"'",$2,$3,".","+",".","ID=InterSpRep_"j ;}' >> $out_gff_cov
## Clean up
rm $temp_bed_sorted
rm $temp_bed
# End of file
