#!/bin/bash
#
# Runs AlphaFold for all monomers .pdb files in a given directory

# Check number of arguments
if [ $# -lt 3 ] || [ $# -gt 3 ] ; then
	echo "Use: AF_mon.sh target_directory output_directory type_of_protein"
	echo "type_of_protein = monomer if you want AlphaFold to run in monomer mode"
	echo "type_of_protein = multimer if you want AlphaFold to run in multimer mode"
	exit
fi

# Check if a correct type of protein is given
if [ $3 != "monomer" ] && [ $3 != "multimer" ] ; then
	echo "Error: no valid type of protein $3"
	echo "Valid types: monomer / multimer"
	exit
fi

# Check if input directory exists
if [ ! -d $1 ] ; then
	echo "Error: directory $1 does not exist"
	exit
fi

# Create output directory
if [ ! -d $2 ] ; then
	echo "Creating directory $2"
	mkdir $2
fi


# Main procedure: run AlphaFold2 for each fasta file.
FILES=$(find $1 -type f -name '*.fasta')

for FILE in $FILES ; do
	FILENAME="${FILE##*/}"
	BASENAME="${FILENAME%.fasta}"
	echo "Processing $FILENAME ..."

	FULLPATH=$(readlink -f "$FILE")
	FULLDIR=$(readlink -f "$2")

	OUTDIR="$2/$BASENAME"

	# Do not run AlphaFold2 if results are already available
	if [ $3 == "monomer" ] ; then

		if [ -d "$OUTDIR" ] && \
		[ -f "$OUTDIR/ranked_0.pdb" ] && \
		[ -f "$OUTDIR/Error_${BASENAME}.txt" ] ; then

			echo "Skipping file $FILENAME (results already available in $OUTDIR)"
			echo -e "\n"
			continue
		fi

	else
		if [ -d "$OUTDIR" ] && \
                [ -f "$OUTDIR/ranked_0.pdb" ] && \
                [ -f "$OUTDIR/Error_${BASENAME}.txt" ] && \
		[ -f "$OUTDIR/${BASENAME}_ipTM.txt" ] ; then

			echo "Skipping file $FILENAME (results already available in $OUTDIR)"
                        echo -e "\n"
                        continue
		fi
	fi

	# Do not run AlphaFold2 if complex is too large (> 4 Monomers)
	if [ $3 == "multimer" ] ; then
		COUNT=$(grep -oE 'PP_[0-9]+' <<< "$FILENAME" | wc -l)
		if [ $COUNT -gt 4 ] ; then
			echo "Skipping $FILENAME. Complex too large ..."
			continue
		fi
	fi

	# Delete all previous files if needed
	if [ -d $OUTDIR ] ; then
		if [ "$(ls -A "$OUTDIR")" ] ; then
			rm -rf "$OUTDIR"/*
			echo "Deleting previous files ..."
		fi
	fi

	# Run AlphaFold2
	mkdir -p "$OUTDIR"
	echo "Running AlphaFold2 for $FILENAME ..."

	python3 ~jnovoa/Software/alphafold/docker/run_docker.py \
		--fasta_paths=$FULLPATH \
		--max_template_date=2025-05-14 \
		--model_preset=$3 \
		-db_preset=full_dbs \
		--data_dir=/data/Alphafold/ \
		--output_dir=$FULLDIR \
		2> "$OUTDIR/Error_${BASENAME}.txt"

	# Retrieve the reesults
	echo "Obtaining best model for $BASENAME ..."
	python3 AF_BestRanked_and_ipTM_BASH.py $OUTDIR $3

	echo "Done"
	echo -e "\n"
done


