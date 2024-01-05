#!/bin/zsh

show_usage() {
    echo "Usage: ./run_ngs.sh [-f <frameshift>] [-n <Path to the reference sequence>] <File with samples>"
    echo "Example: ./run_ngs.sh -f 1 -n ref_seq.fasta samples.txt"
}

error_exit() {
    echo "Error: $1" >&2
    exit "${2:-1}"
}

while getopts ":f:n:" opt; do
    case $opt in
        f) frameshift="$OPTARG";;
        n) path_to_nat_sequence="$OPTARG";;
        \?) echo "Invalid argument: -$OPTARG" >&2; show_usage; exit 1;;
        :) echo "argument -$OPTARG expects defintion." >&2; show_usage; exit 1;;
    esac
done

# Move the position of the argument to take the "optind" into account -> $1 = samples
shift $((OPTIND - 1))
sample_file="$1"

# Checks whether -f and -n are defined
if [ -z "$frameshift" ] || [ -z "$path_to_nat_sequence" ]; then
    show_usage
    error_exit "The options -n and -f must be defined." 1
fi

# Checks whether the number of arguments is correct
if [ "$#" -ne 1 ]; then
    show_usage
    error_exit "Incorrect number of arguments. Please enter frameshift, path to the native sequence and sample file." 1
fi

mkdir -p ./demultiplex || error_exit "Demultiplex directory could not be created."

# Demultiplex and .seq to .fasta
echo "Starte Demultiplexing..."
while IFS=',' read -r fasta sample bc_start bc_end length patch || [[ -n "$fasta" ]]; do
    if [[ -n "$fasta" ]]; then
        # Demultiplex
        grep -o "$bc_start.*$bc_end" "$fasta.fasta" | sed -e "s/.*$bc_start\(.*\)$bc_end.*/\1/" | awk "length(\$0) ==$length" > "demultiplex/${sample}.seq" || error_exit "Demultiplex-Fehler für $sample."
        
        # .seq to .fasta
        SeqFile="demultiplex/${sample}.fasta"

        if [[ -e "$SeqFile" ]]; then
            error_exit "The file '$SeqFile' already exists." 2
        else
            counter=0
            while IFS= read -r line; do
                let counter++
                echo ">${sample}_${counter}" >> "$SeqFile"
                echo "$line" >> "$SeqFile" 
            done < "demultiplex/${sample}.seq"
            seq_count=$(wc -l < "$SeqFile")
            echo "${sample} demultiplexed - Number of sequnces: ${seq_count}"
        fi
    fi
done < "$sample_file"

# Ergebnisse erstellen
echo "Starte Plotting..."
mkdir -p ./results || error_exit "Results directory could not be created."
python plot.py "$frameshift" "$path_to_nat_sequence" "$sample_file" || error_exit "Python script error."
